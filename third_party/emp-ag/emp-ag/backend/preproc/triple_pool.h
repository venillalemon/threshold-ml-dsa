#ifndef EMP_AG_BACKEND_PREPROC_TRIPLE_POOL_H__
#define EMP_AG_BACKEND_PREPROC_TRIPLE_POOL_H__
#include "emp-ag/backend/preproc/auth_share_pool.h"
#include "emp-ag/backend/crypto/fs_hash.h"  // FsHash (FS/check digests)
#include "emp-ag/backend/crypto/feq.h"     // fs_coin_pair, symmetric F_eq
#include "emp-ag/backend/crypto/fzero.h"   // fzero_xor
#include "emp-ag/backend/crypto/protocols.h"   // sampleRandom
#include "emp-ag/backend/netmp.h"
#include "emp-ag/backend/profiling.h"
#include "emp-ag/backend/runtime.h"        // parallel_for, lane_range
#include "emp-ag/backend/crypto/bucketing.h"   // emp::ag::bucket_size (shared)
#include "emp-ag/backend/crypto/auth_open.h"   // local authenticated open checks
#include "emp-tool/runtime/crypto/mitccrh.h"
#include <climits>
#include <cstdint>
#include <thread>
#include "emp-tool/runtime/core/block.h"
#include "emp-tool/runtime/crypto/prg.h"      // PRG
#include <vector>
#include "emp-tool/runtime/crypto/hash.h"
#include <algorithm>
#include <cstring>
#include <functional>
#include <memory>

namespace emp::ag {

// ---- Bucketing leaves used only by TriplePool: the cyclic-shift fold, the batched
//      schedule (run_bucketing_core), and the windowed F_eq. bucket_size stays in
//      backend/crypto/bucketing.h because its pure-math unit test (test/bucketing_test.cpp)
//      includes it standalone.

// ===== cyclic-shift bucket fold =====
// Source index of layer element i under cyclic shift r in [0,L): the wraparound
// splits at i==L-r. THE byte-identity tripwire — an off-by-one here drifts every
// downstream transcript — so it is single-sourced and used by all three leaves.
inline int bucket_fold_src(int i, int L, int r) {
  const int cut = L - r;
  return (i < cut) ? i + r : i + r - L;
}

// Draw an exactly uniform cyclic-shift index. A plain `word % L` is biased
// unless L divides the word range; widening that word to 128 bits only makes
// the mismatch tiny, while the bucketing proof assumes a uniform rotation.
// Lemire-style threshold rejection over a fixed little-endian uint64_t is exact,
// deterministic across architectures, and rejects with probability < 2^-33 for
// the supported L<=INT_MAX. Each attempt consumes one PRG block, just like the
// old four-byte draw's small-buffer path.
inline int bucket_uniform_shift(PRG &prg, int L) {
  expecting(L > 0, "bucket_uniform_shift: L must be positive");
  const uint64_t bound = (uint64_t)L;
  const uint64_t threshold = (uint64_t(0) - bound) % bound;
  uint64_t raw;
  do {
    unsigned char bytes[8];
    prg.random_data_unaligned(bytes, sizeof(bytes));
    raw = 0;
    for (int i = 0; i < 8; ++i) raw |= (uint64_t)bytes[i] << (8 * i);
  } while (raw < threshold);
  return (int)(raw % bound);
}

// Phase 1: a/c XOR fold. acc.a ^= layer.a[src]; acc.c ^= layer.c[src] over the
// (mac,key) pair; b is untouched. No d, no IO. [i0, i1) bounds the OUTPUT index
// range so the caller can width-split one peer's fold across pool workers: the
// per-element XOR-RMW makes ranges of the same (peer, layer) race-free, while
// two LAYERS of the same peer must stay ordered (both RMW the same acc slots).
inline void bucket_fold_pre_ac(block *am, block *ak, const block *lm,
                               const block *lk, int L, int r, int i0, int i1) {
  for (int i = i0; i < i1; ++i) {
    const int src = bucket_fold_src(i, L, r);
    am[i]         = am[i]         ^ lm[src];
    am[2 * L + i] = am[2 * L + i] ^ lm[2 * L + src];
    ak[i]         = ak[i]         ^ lk[src];
    ak[2 * L + i] = ak[2 * L + i] ^ lk[2 * L + src];
  }
}

// The public d-bit: bit0 of the b-region MAC of acc vs the shifted layer. Filled
// ONCE from the chosen peer slot (b is never written by the a/c fold, so it may
// run concurrently with any bucket_fold_pre_ac range, including that peer's own).
// The caller exchanges d over its own transport and combines (XOR for 2PC; the
// cross-peer fold for n-party). [i0, i1) as above.
inline void bucket_fold_dbit(const block *am, const block *lm, int L, int r,
                             unsigned char *d_out, int i0, int i1) {
  for (int i = i0; i < i1; ++i) {
    const int src = bucket_fold_src(i, L, r);
    d_out[i] = (unsigned char)(getLSB(am[L + i]) ^ getLSB(lm[L + src]));
  }
}

// Phase 3: masked c-correction, keyed on the ALREADY-COMBINED public d (d[1] after the
// cross-peer fold; d_me^d_peer at nP==2). d_comb is read-only and shared across peers.
// [i0, i1) as above.
inline void bucket_fold_post(block *am, block *ak, const block *lm,
                             const block *lk, int L, int r,
                             const unsigned char *d_comb, int i0, int i1) {
  for (int i = i0; i < i1; ++i) {
    const int src = bucket_fold_src(i, L, r);
    const block m = select_mask[d_comb[i]];
    am[2 * L + i] = am[2 * L + i] ^ (lm[src] & m);
    ak[2 * L + i] = ak[2 * L + i] ^ (lk[src] & m);
  }
}

// ===== batched bucketing schedule =====
template <class BucketPolicy>
void run_bucketing_core(BucketPolicy &pol, int B, int L) {
  // ---- Pass 1: generate + COMMIT all B layers. nP==2 batches their row flight
  // and dependent s flight; nP>=3 retains its existing per-layer construction.
  // No coin yet, so the entire B*L set is fixed before draw_coins. ----
  pol.generate_layers();

  // ---- The ONE post-commitment coin -> the B-1 cyclic shifts. One PRG stream from
  // one seed (so the shifts are jointly pseudorandom). The PRG outlives the shift
  // draw and is handed to finalize_check: nP>=3 CONTINUES this exact stream for its
  // batched-RLC challenge phi_cat (so the challenge is drawn from the same post-commitment
  // coin, right after the shifts); nP==2 ignores it (feq_check). ----
  block coin = pol.draw_coins();
  PRG prg(&coin);
  std::vector<int> shift(B, 0);
  for (int k = 1; k < B; ++k) shift[k] = bucket_uniform_shift(prg, L);

  // ---- Pass 2: shifts known -> fold each layer's a/c + the b-region d-bit (local),
  // exchange ALL (B-1)*L d-bits in ONE round, then apply the masked c-correction. ----
  for (int k = 1; k < B; ++k) pol.fold_pre_and_dbit(k, shift[k]);
  pol.exchange_all_d();
  for (int k = 1; k < B; ++k) pol.fold_post(k, shift[k]);

  // ---- ONE terminal check over all B layers (the batched F_eq / RLC zero-sum).
  pol.finalize_check(prg);
}

// ===== windowed F_eq (leaky-AND) =====
// Number of low bytes of each half-gate row actually transmitted/checked.
inline int laand_window_bytes(int ssp) {
  return std::min((int)sizeof(block), (ssp + 2 + 7) / 8);
}

// Block whose low laand_window_bytes(ssp) bytes are 0xFF and the rest 0 — used by
// the nP>=3 projected alpha check so both ends carryless-multiply the same
// windowed object (the high bytes are not reconciled cross-party once G is
// truncated; both ends zero them).
inline block laand_window_mask(int ssp) {
  const int gb = laand_window_bytes(ssp);
  unsigned char mb[sizeof(block)] = {0};
  for (int i = 0; i < gb; ++i) mb[i] = (unsigned char)0xFF;
  block m;
  memcpy(&m, mb, sizeof(block));
  return m;
}

// Fixed stripe width (in elements) for the parallel F_eq window fold below. It
// shapes the CHECK VALUE (stripe digests are what the running hash binds), so it
// must be identical at both parties: a compile-time constant, never a function of
// the local pool size.
inline constexpr int kFeqStripeElems = 1 << 16;

// Fold one bucket layer's L-vector into the running F_eq hash, restricted to the
// low `gb` bytes of each element (the transmitted window). The fold is a fixed
// TWO-LEVEL hash: the window bytes are gathered + hashed per kFeqStripeElems-element
// stripe — independent stripes, so the gather and the hashing run pool-parallel;
// the flat single-stream fold was the last big serial block on the nP==2
// draw+seed critical path — and the ordered stripe digests feed the running feq.
// This binds exactly the same window bytes (collision resistance composes across
// the levels); the check VALUE differs from the flat fold, so both parties carry
// the same stripe constant and derive the same stripe count from the common L.
// One terminal symmetric F_eq fires per bucket; the caller owns it.
inline void accumulate_leaky_and_window(FsHash &feq, const block *src, int L, int gb,
                                        std::vector<unsigned char> &scratch,
                                        ThreadPool *pool) {
  // Caller-owned grow-only scratch (a plain member, passed by ref -- no per-access
  // overhead): this runs B x per chunk, so a fresh L*gb buffer each
  // call is mmap/fault/munmap churn. Stripes write disjoint ranges.
  const size_t need = (size_t)L * gb;
  if (scratch.size() < need) scratch.resize(need);
  const int ns = (L + kFeqStripeElems - 1) / kFeqStripeElems;
  std::vector<char> dg((size_t)ns * FsHash::DIGEST_SIZE);
  vector<future<void>> res;
  for (int s = 0; s < ns; ++s) {
    res.push_back(pool->enqueue([&, s]() {
      const int k0 = s * kFeqStripeElems;
      const int k1 = std::min(L, k0 + kFeqStripeElems);
      unsigned char *dst = scratch.data() + (size_t)k0 * gb;
      for (int k = k0; k < k1; ++k)
        memcpy(scratch.data() + (size_t)k * gb, &src[k], gb);
      FsHash::hash_once(dg.data() + (size_t)s * FsHash::DIGEST_SIZE, dst,
                      (size_t)(k1 - k0) * gb);
    }));
  }
  joinNclean(res);
  feq.put(dg.data(), dg.size());
}

// Batched nP==2 production already keeps each low-byte window contiguous, so
// preserve the exact same stripe-hash construction without gathering blocks
// through a second scratch buffer.
inline void accumulate_leaky_and_window_bytes(FsHash &feq,
                                              const unsigned char *src,
                                              int L, int gb, ThreadPool *pool) {
  const int ns = (L + kFeqStripeElems - 1) / kFeqStripeElems;
  std::vector<char> dg((size_t)ns * FsHash::DIGEST_SIZE);
  vector<future<void>> res;
  for (int s = 0; s < ns; ++s) {
    res.push_back(pool->enqueue([&, s]() {
      const int k0 = s * kFeqStripeElems;
      const int k1 = std::min(L, k0 + kFeqStripeElems);
      FsHash::hash_once(dg.data() + (size_t)s * FsHash::DIGEST_SIZE,
                        src + (size_t)k0 * gb,
                        (size_t)(k1 - k0) * gb);
    }));
  }
  joinNclean(res);
  feq.put(dg.data(), dg.size());
}


// Π_LaAND followed by Π_Prep bucketing + amortized pool: AND-triple
// generation via half-gate garbled-AND + phi-based triple check (P:LaAND),
// then circular-shift bucketing to remove leakage (P:aAND).
//
// Deviations from the paper:
//
// (1) FZero (P:LaAND step 3) is implemented via pairwise seed-and-expand
//     under the ICM assumption on the AES-based PRG: each pair shares a
//     λ-bit seed, both expand into LB+1 blocks via PRG, XOR into local
//     (z^i, u^i). Each seed contributes to exactly two parties so
//     Σ_p z^p_k = 0 and Σ_p u^p = 0. (Same construction as
//     AuthSharePool::check1.)
//
// (2) F_Com is collapsed to a single send for both d^i (steps 6→7) and
//     α^i (steps 9→10).  No value derived from the committed bits is
//     sent before the matching open, so equivocation is impossible —
//     the commitment phase is redundant.  α^i is hash-committed first
//     (one round) and then opened, since the α-check needs binding.
//
// (3) The B-1 sequential OpenToAll calls in P:aAND step 4 (one per
//     non-first leaky triple in each bucket) are batched into a single
//     all-to-all exchange of the full d-vector, dropping the bucket
//     combine to one round.
//
// This is the per-call triple / AND-share generator (no pool): compute_inplace
// produces the AND gates' σ = λ_α∧λ_β shares for one circuit pass via the
// layered leaky-AND + bucketing (one in-place representative layer + B-1
// out-of-place sacrifice layers).
template <int nP>
class TriplePool {
public:
  ThreadPool *pool;
  int party;
  int ssp;              // statistical security parameter (default 40)
  NetIOMP<nP> *io;
  block session_id;
  AuthSharePool<nP> abit;
  block Delta;
  uint64_t checked_open_epoch_ = 0;
  BlockVec originalMAC_[nP + 1], originalKEY_[nP + 1];
#ifdef EMP_AG_LOCAL_TEST_HOOKS
  // Unit-test-only malformed-opening injection; absent from production.
  std::function<void(const char *, unsigned char *, size_t)> local_open_test_hook;
#endif
  PRG prg;
  // Half-gate hash seed (MITCCRH start_point), piggybacked off the per-chunk
  // post-commitment bucket coin in BucketPolicy::draw_coins (one extra domain-separated
  // squeeze, no dedicated coin round). AGMPCProtocol::draw_and_seed reads it to seed
  // ctx.mitc after compute_inplace; refreshed on every compute_inplace that has ANDs.
  block halfgate_seed_ = zero_block;
  // Batched LaAND alpha-check stash (nP>=3): each layer's windowed RLC vector +
  // its FZero u_zero are accumulated here across all B layers of a bucket; the ONE
  // B*L RLC zero-sum runs in compute_inplace (INC-2). Cleared per compute_inplace.
  BlockVec alpha_masked_cat_;
  std::vector<block> alpha_uzero_;
  // Reused per-chunk bucketing scratch (grow-only). compute_inplace runs once per
  // protocol chunk; these per-peer
  // accumulator + rep-COT buffers and the B-1 sacrifice layers (3*L each) + d-vectors
  // are kept resident across chunks instead of being mmap/first-touch-fault/munmap'd
  // every chunk. Bytes overwritten each chunk; grow-only => no realloc once the
  // steady-state chunk size is reached.
  BlockVec bkt_accMAC_[nP + 1], bkt_accKEY_[nP + 1];
  BlockVec bkt_rMAC_[nP + 1], bkt_rKEY_[nP + 1];
  std::vector<std::vector<BlockVec>> bkt_sac_macs_, bkt_sac_keys_;
  std::vector<unsigned char> bkt_bd_[nP + 1];
  // produce_leaky_ands_halfgate scratch (per-layer, B x per chunk; grow-only). All
  // assigned/explicitly-zeroed before read except z_u (XOR-accumulated by fzero_xor
  // -> caller re-zeros via std::fill each call). + the per-chunk RLC challenge phi_cat
  // (nP>=3 finalize_check) and the compute_inplace reveal buffers x/y.
  BlockVec bkt_tKEYphi_[nP + 1], bkt_tMACphi_[nP + 1];
  BlockVec bkt_phi_, bkt_z_u_, bkt_phi_cat_;
  std::vector<unsigned char> bkt_s_[nP + 1];
  std::vector<unsigned char> bkt_x_[nP + 1], bkt_y_[nP + 1];
  std::vector<unsigned char> bkt_window_;   // accumulate_leaky_and_window scratch
  BlockVec bkt_fzero_;                       // fzero_xor contribution scratch
  // grp (optional): a NetIOMPGroup (mesh(0) == io) enabling lane-parallel COT in
  // AuthSharePool AND the leaky-AND garbled-AND row transfer here; null => single
  // mesh. L_lanes = the shared lane count.
  NetIOMPGroup<nP> *grp = nullptr;
  int L_lanes = 1;
  TriplePool(NetIOMP<nP> *io, ThreadPool *pool, int party, int ssp = 40,
             NetIOMPGroup<nP> *grp_in = nullptr,
             block session_id_in = zero_block,
             BaseOtKind base_ot_kind = BaseOtKind::CSW,
             int64_t planned_sender_cots = -1,
             int64_t planned_receiver_cots = -1,
             const block *supplied_delta = nullptr)
      : pool(pool), party(party), ssp(ag::checked_statistical_security(ssp)), io(io),
        session_id(session_id_in),
        abit(io, pool, party, nullptr, grp_in, session_id_in, base_ot_kind,
             planned_sender_cots, planned_receiver_cots, supplied_delta),
        Delta(abit.Delta),
        grp(grp_in), L_lanes(grp_in ? std::max(1, grp_in->size()) : 1) {}

  // The mesh carrying lane t (lane 0 == io == grp->mesh(0)); see AuthSharePool.
  NetIOMP<nP> *lane_mesh(int t) { return (t == 0) ? io : &grp->mesh(t); }

  // Canonical per-(peer,lane) transfer fan-out over the group meshes, into ONE join
  // barrier. Splits [0,n) into T=min(L_lanes,n) contiguous ag::lane_range slices; lane t
  // rides mesh t (lane 0 == mesh 0). `send_slice(ch,peer,lo,hi)` and `recv_slice(...)` do
  // the transfer for that slice on channel ch. Order = peer-major, lanes 0..T-1, send+recv
  // adjacent — THE deadlock-safe enqueue contract, defined once here for every bulk open
  // (d-vector, s-open, x/y). n==0 is a no-op; T==1 keeps everything on mesh 0 (byte-identical).
  template <class SendF, class RecvF>
  void lane_fan_out(int n, SendF &&send_slice, RecvF &&recv_slice) {
    const int T = std::min(L_lanes, std::max(1, n));
    vector<future<void>> res;
    for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
      for (int t = 0; t < T; ++t) {
        auto pr = ag::lane_range(t, n, T);
        NetIOMP<nP> *ch = lane_mesh(t);
        res.push_back(pool->enqueue([&send_slice, ch, peer, pr]() { send_slice(ch, peer, pr.first, pr.second); }));
        res.push_back(pool->enqueue([&recv_slice, ch, peer, pr]() { recv_slice(ch, peer, pr.first, pr.second); }));
      }
    }
    joinNclean(res);
  }

  // The two NetIO channels behind NetIOMP for one peer (see netmp.h: the
  // smaller party drives `ios`, the larger drives `ios2`). send_data/recv_data
  // pick these; we reach them directly to send bit vectors via the channel's
  // packed send_bool/recv_bool (8 bits/byte) instead of one byte per bit.
  NetIO *send_ch(int dst) { return party < dst ? io->ios[dst] : io->ios2[dst]; }
  NetIO *recv_ch(int src) { return src < party ? io->ios[src] : io->ios2[src]; }

  // Bucket size B vs. number of triples ℓ — derived from the proven aAND residual-leakage
  // bound (2B+1)/ℓ^{(B-1)} <= 2^{-ssp} (the aAND theorem), not tuned by hand (a wrong B
  // silently widens residual leakage, invisible to every functional test). The same bound
  // governs every party count. See emp-ag/backend/crypto/bucketing.h for the derivation.
  int get_bucket_size(int size) {
    return ag::bucket_size(size, ssp);
  }
  // Multi-party half-gate leaky-AND over one layer of L candidates. The layer is
  // laid out as a=[0,L), b=[L,2L), r/c=[2L,3L). On return r has been folded into
  // c=a AND b in place, and the LaAND alpha check has been run for this layer.
  void produce_leaky_ands_halfgate(BlockVec tMAC[nP + 1],
                                   BlockVec tKEY[nP + 1], int L, FsHash &feq) {
    // ONE leaky-AND generator for every party count. The nP==2 collapse is a set of
    // conditioned branches, not a separate construction: FZero + the cross-party phi
    // sum are nP>2 only (at two parties the blind C is local), and the half-gate row
    // uses MITCCRH hash_cir at nP==2 / hash at nP>2. Both fold their per-peer rows into
    // tKEYphi[party] (the shared result buffer), so the s-open, the d-combine, and the
    // terminal-check accumulation below are byte-identical across party counts.
    BlockVec(&tKEYphi)[nP + 1] = bkt_tKEYphi_;
    std::vector<unsigned char>(&s)[nP + 1] = bkt_s_;
    const int any_peer = (party == 1) ? 2 : 1;
    const int gb = ag::laand_window_bytes(ssp);
    block u_zero = zero_block;   // FZero blind (nP>2 only; bound by the batched check)
    for (int i = 0; i <= nP; ++i)
      if (s[i].size() < (size_t)L) s[i].resize((size_t)L);
    if (tKEYphi[party].size() < (size_t)L) tKEYphi[party].resize((size_t)L);

#ifdef AG_PROFILE
    uint64_t _tp_hg0 = agmpc_now_ns();
#endif
    // --- phi: the half-gate row correction (the blind C, generalized). At nP>2 it
    // starts from the FZero zero-share (carrying u_zero, bound by the batched check);
    // at nP==2 there is no FZero (z_u=0), so phi reduces to the bare local blind C. Then
    // add the (sel&Delta) term + the cross-peer key/mac sum (one peer at nP==2). ---
    BlockVec &phi = bkt_phi_;
    if (phi.size() < (size_t)L) phi.resize((size_t)L);
    if constexpr (nP > 2) {
      BlockVec &z_u = bkt_z_u_;
      if (z_u.size() < (size_t)(L + 1)) z_u.resize((size_t)(L + 1));
      std::fill(z_u.begin(), z_u.begin() + (L + 1), zero_block);  // fzero_xor XOR-accumulates
      fzero_xor<nP>(io, &prg, pool, party, z_u.data(), L + 1, bkt_fzero_);
      for (int k = 0; k < L; ++k) phi[k] = z_u[k];
      u_zero = z_u[L];
    } else {
      std::fill(phi.begin(), phi.begin() + L, zero_block);
    }
    int width = (L + pool->size() - 1) / pool->size();
    {
      vector<future<void>> res;
      for (int wi = 0; wi < (int)pool->size(); ++wi) {
        int st = wi * width, ed = std::min((wi + 1) * width, L);
        res.push_back(pool->enqueue([this, st, ed, L, any_peer, &tKEY, &tMAC, &phi]() {
          for (int k = st; k < ed; ++k) {
            phi[k] ^= (select_mask[getLSB(tMAC[any_peer][L + k])] & Delta);
            for (int j = 1; j <= nP; ++j) if (j != party) {
              phi[k] ^= tKEY[j][L + k];
              phi[k] ^= tMAC[j][L + k];
            }
          }
        }));
      }
      joinNclean(res);
    }

    // --- per-peer fused half-gate (MITCCRH hash_cir under the per-pair FS seed): one
    // peer at nP==2, n-1 in parallel at nP>2. The garbler ships row = H(key)^H(key^Delta)^phi
    // (low gb bytes; F_eq / the windowed alpha check covers exactly that window). The
    // evaluator hashes [mac,key] -> HM,HK and folds the per-peer term HK ^ (HM ^ (w&sel))
    // ^ key[2L] ^ mac[2L] into tMACphi[peer]. ---
    BlockVec(&tMACphi)[nP + 1] = bkt_tMACphi_;
    for (int i = 1; i <= nP; ++i)
      if (i != party && tMACphi[i].size() < (size_t)L) tMACphi[i].resize((size_t)L);
    const int hg_chunk = 1 << 16;
#ifdef AG_PROFILE
    int64_t _phi0 = io->count();
#endif
    {
      // Per-lane CCR key salt in the FREE low 64:
      // MITCCRH derives its AES key as start_point ^ makeBlock(gid>>ReuseShift, 0) -- the gid
      // COUNTER marches through the HIGH 64 bits of the key, the LOW 64 are always 0 (free).
      // So salt each lane's seed in that free low 64: setS(seed ^ makeBlock(0, lane)) => key =
      // seed ^ makeBlock(bucket, lane) -- counter in the high 64, lane in the low 64, a
      // DISJOINT field, unique per (lane, bucket), no key reuse. Blackbox, one call, no
      // gid-range constraint; lane 0 salts by 0 => byte-identical to the pre-lane build at
      // T==1. The AND-index range [0,L) is split into T=min(L_lanes,L) contiguous slices; lane
      // t transfers slice t on mesh t (lane 0 == mesh 0), assembling into the disjoint tMACphi
      // range so the check folds on the full [0,L) unchanged. Only the bulk row bytes move;
      // control-message flights stay on mesh 0, while local FS challenges snapshot every
      // lane transcript before deriving one shared seed (R6/R7).

      // Sender half-gate row transfer for AND-index range [lo,hi) under lane `lane`, on ch.
      auto laand_send = [this, &tKEY, &phi, hg_chunk, gb](NetIOMP<nP> *ch, int peer, block seed, int lo, int hi, int lane) {
        MITCCRH<8> mitc;
        mitc.setS(seed ^ makeBlock(0, (uint64_t)lane));   // lane in the free low 64 of the key
        block pad[16];
        const int len = hi - lo;
        std::vector<unsigned char> gbuf((size_t)std::min(hg_chunk, len) * gb);
        for (int ci = 0; ci < (len + hg_chunk - 1) / hg_chunk; ++ci) {
          int st = lo + hg_chunk * ci, ed = std::min(lo + hg_chunk * (ci + 1), hi);
          for (int k0 = st; k0 < ed; k0 += 8) {
            int batch = std::min(8, ed - k0);
            for (int j = 0; j < 8; ++j) {
              block kk = (j < batch) ? tKEY[peer][k0 + j] : zero_block;
              pad[2 * j]     = kk;
              pad[2 * j + 1] = kk ^ Delta;
            }
            mitc.hash_cir<8, 2>(pad);
            for (int j = 0; j < batch; ++j) {
              block row = pad[2 * j] ^ pad[2 * j + 1] ^ phi[k0 + j];   // full G row
              memcpy(gbuf.data() + (size_t)(k0 + j - st) * gb, &row, gb); // ship low gb bytes
            }
          }
#ifdef EMP_AG_TEST_FLIP_LAAND
          // TEST-ONLY (never in production): a malicious party corrupts ONE shipped leaky-AND
          // row (party 2's row to party 1), so P1's tKEYphi diverges. The leaky-AND consistency
          // check -- F_eq at nP==2, the batched projected-alpha RLC at nP>=3 -- must abort during
          // preprocessing. It must be a SINGLE row: flipping symmetrically at every sender would
          // cancel in the XOR-combined tKEYphi at nP>=3 (each receiver gets the same delta from
          // all peers), a false-negative that hides the fault.
          if (ci == 0 && ed > st && party == 2 && peer == 1 && lane == 0) gbuf[0] ^= 0x1;
#endif
          ch->send_data(peer, gbuf.data(), (size_t)(ed - st) * gb);
          ch->flush(peer);
        }
      };

      // Receiver: recompute tMACphi[peer][k] for k in [lo,hi) from the shipped rows on ch.
      auto laand_recv = [&tMAC, &tKEY, &tMACphi, any_peer, L, hg_chunk, gb](NetIOMP<nP> *ch, int peer, block seed, int lo, int hi, int lane) {
        MITCCRH<8> mitc;
        mitc.setS(seed ^ makeBlock(0, (uint64_t)lane));   // lane in the free low 64 of the key
        block pad[16];
        const int len = hi - lo;
        BlockVec wire(std::min(hg_chunk, len));
        std::vector<unsigned char> gbuf((size_t)std::min(hg_chunk, len) * gb);
        for (int ci = 0; ci < (len + hg_chunk - 1) / hg_chunk; ++ci) {
          int st = lo + hg_chunk * ci, ed = std::min(lo + hg_chunk * (ci + 1), hi);
          ch->recv_data(peer, gbuf.data(), (size_t)(ed - st) * gb);
          for (int r = 0; r < ed - st; ++r) {   // zero-extend each gb-byte row to a block
            wire[r] = zero_block;
            memcpy(&wire[r], gbuf.data() + (size_t)r * gb, gb);
          }
          for (int k0 = st; k0 < ed; k0 += 8) {
            int batch = std::min(8, ed - k0);
            for (int j = 0; j < 8; ++j) {
              pad[2 * j]     = (j < batch) ? tMAC[peer][k0 + j] : zero_block;  // -> H(M)
              pad[2 * j + 1] = (j < batch) ? tKEY[peer][k0 + j] : zero_block;  // -> H(K)
            }
            mitc.hash_cir<8, 2>(pad);
            for (int j = 0; j < batch; ++j) {
              int k = k0 + j;
              block HM = pad[2 * j], HK = pad[2 * j + 1];
              block E = HM ^ (wire[k - st] & select_mask[getLSB(tMAC[any_peer][k])]);
              tMACphi[peer][k] = HK ^ E ^ tKEY[peer][2 * L + k] ^ tMAC[peer][2 * L + k];
            }
          }
        }
      };

      // The leaky-AND transfer keeps its own fan-out (one per-peer seed shared by all lanes,
      // and the measurement cap below) rather than lane_fan_out, but uses the shared
      // ag::lane_range.
      int T = std::min(L_lanes, std::max(1, L));   // T_eff = min(lanes, L)  (R8: no empty slice)
#ifdef EMP_AG_LAAND_LANE_CAP
      // Measurement-only: cap the leaky-AND lane count independently of the COT lanes
      // (build -DEMP_AG_LAAND_LANE_CAP=1 to keep COT lane-split but run the leaky-AND
      // transfer single-mesh, isolating this phase's lane contribution). No default effect.
      T = std::min(T, EMP_AG_LAAND_LANE_CAP);
#endif

      vector<future<void>> res;
      for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
        // Snapshot both directed transcripts from every mesh, then share the resulting
        // per-peer seed across lanes (the lane id further domain-separates the CCR key).
        block seed = fs_coin_pair_group<nP>(
            io, grp, peer, "AGMPC laand seed", session_id);
        for (int t = 0; t < T; ++t) {   // canonical order: peer-major, lanes 0..T-1, dir-adjacent (R4)
          auto [lo, hi] = ag::lane_range(t, L, T);
          NetIOMP<nP> *ch = lane_mesh(t);
          res.push_back(pool->enqueue([&laand_send, ch, peer, seed, lo, hi, t]() { laand_send(ch, peer, seed, lo, hi, t); }));
          res.push_back(pool->enqueue([&laand_recv, ch, peer, seed, lo, hi, t]() { laand_recv(ch, peer, seed, lo, hi, t); }));
        }
      }
      joinNclean(res);   // all-lane flush barrier: every lane's send+recv joined before draw_coins (R5)
    }
#ifdef AG_PROFILE
    g_ag_phi_bytes += (uint64_t)(io->count() - _phi0);
#endif

    // --- combine: sum the per-peer terms (one peer at nP==2) + the global phi/Delta
    // correction -> the leaky-AND result in tKEYphi[party]. ---
    {
      vector<future<void>> res;
      for (int wi = 0; wi < (int)pool->size(); ++wi) {
        int st = wi * width, ed = std::min((wi + 1) * width, L);
        res.push_back(pool->enqueue([this, st, ed, L, any_peer, &tKEYphi, &tMACphi, &tMAC, &phi]() {
          for (int k = st; k < ed; ++k) {
            block t = zero_block;
            for (int j = 1; j <= nP; ++j) if (j != party) t = t ^ tMACphi[j][k];
            t = t ^ (phi[k] & select_mask[getLSB(tMAC[any_peer][k])]);
            t = t ^ (Delta & select_mask[getLSB(tMAC[any_peer][2 * L + k])]);
            tKEYphi[party][k] = t;
          }
        }));
      }
      joinNclean(res);
    }
#ifdef AG_PROFILE
    g_ag_tp_hg_ns += agmpc_now_ns() - _tp_hg0;
    uint64_t _tp_so0 = agmpc_now_ns();
#endif

    // s-open (shared): d = LSB1(result). Exchange this party's s with every peer.
    // Lane-split (§A.3, LANE_CLEAN): slice [0,L) across the meshes. Pre-draw_coins, but
    // covered by the SAME all-bytes-on-wire-before-draw_coins ordering as the leaky-AND
    // rows above (this join completes before produce_leaky_ands_halfgate returns) -- the
    // SEAM-B post-commitment resolution applies identically. Combine stays on assembled s[].
    {
      for (int k = 0; k < L; ++k)
        s[party][k] = LSB1(tKEYphi[party][k]);
      lane_fan_out(L,
        [this, &s](NetIOMP<nP> *ch, int peer, int lo, int hi) {
          ch->send_bool(peer, (const bool *)(s[party].data() + lo), hi - lo);
          ch->flush(peer);
        },
        [&s](NetIOMP<nP> *ch, int peer, int lo, int hi) {
          ch->recv_bool(peer, (bool *)(s[peer].data() + lo), hi - lo);
        });
    }

#ifdef AG_PROFILE
    g_ag_tp_sopen_ns += agmpc_now_ns() - _tp_so0;
    uint64_t _tp_ac0 = agmpc_now_ns();
#endif
    // d-combine (shared): fold d = ⊕_p s_p into region 2 (r -> c) and the result
    // vector. Width-split: every write is k-indexed (s[0][k], the 2L+k share slots,
    // T[k]) and every read is done, so ranges are race-free and the split is
    // byte-identical — this loop was a serial block on the draw+seed critical path.
    {
      block *T = tKEYphi[party].data();
      block dxor = Delta ^ bit0_mask;
      ag::parallel_for(pool, L, [&](int lo, int hi) {
        for (int k = lo; k < hi; ++k) {
          s[0][k] = 0;
          for (int i = 1; i <= nP; ++i)
            s[0][k] = (s[0][k] != s[i][k]);
          block mask_s = select_mask[s[0][k]];
          if (party == 1) {
            if (s[0][k]) {
              for (int j = 1; j <= nP; ++j) if (j != party)
                tMAC[j][2 * L + k] = tMAC[j][2 * L + k] ^ bit0_mask;
            }
          } else {
            tKEY[1][2 * L + k] = tKEY[1][2 * L + k] ^ (dxor & mask_s);
          }
          T[k] = T[k] ^ (Delta & mask_s);
        }
      });
    }

    // Terminal-check accumulation: batched per bucket. Each layer accumulates its
    // windowed check object; ONE terminal check fires per bucket in compute_inplace
    // (the Func[Rand]-after-Gen ordering -- soundness Delta 1/2).
    if constexpr (nP == 2) {
      // F_eq: fold the low-gb window of this layer's L-vector into the running hash;
      // one terminal feq_check fires in compute_inplace. The hash over the whole
      // windowed vector IS the check.
      ag::accumulate_leaky_and_window(feq, &tKEYphi[party][0], L, gb, bkt_window_, pool);
    } else {
      // nP>=3: STASH this layer's windowed RLC vector (the high bytes are not
      // reconciled cross-party once G is truncated; both ends zero them) + its
      // FZero u_zero. The single B*L GF(2^128) RLC zero-sum runs once per bucket in
      // compute_inplace over ALL B layers -- tighter than B per-layer checks by a
      // factor B, and the challenge is drawn post-commitment (after the whole B*L
      // set is on the wire).
      const block window_mask = ag::laand_window_mask(ssp);
      size_t off = alpha_masked_cat_.size();
      alpha_masked_cat_.resize(off + (size_t)L);
      for (int k = 0; k < L; ++k)
        alpha_masked_cat_[off + k] = tKEYphi[party][k] & window_mask;
      alpha_uzero_.push_back(u_zero);
    }
#ifdef AG_PROFILE
    g_ag_tp_acheck_ns += agmpc_now_ns() - _tp_ac0;
#endif
  }

  // Two-party batched LaAND schedule: every layer's COTs are fixed first, one
  // transcript root domain-separates the B row seeds, all rows cross in one
  // flight, and all dependent s bits cross in the next.  The authenticated
  // layer representation and final F_eq input are unchanged.
  void produce_leaky_ands_halfgate_batch_2pc(
      BlockVec accMAC[nP + 1], BlockVec accKEY[nP + 1],
      std::vector<std::vector<BlockVec>> &sac_macs,
      std::vector<std::vector<BlockVec>> &sac_keys,
      int B, int L, FsHash &feq) {
    if constexpr (nP == 2) {
      const int peer = (party == 1) ? 2 : 1;
      const int gb = ag::laand_window_bytes(ssp);
      const int BL = B * L;
      auto layer_mac = [&](int k) -> BlockVec * {
        return k == 0 ? accMAC : sac_macs[k].data();
      };
      auto layer_key = [&](int k) -> BlockVec * {
        return k == 0 ? accKEY : sac_keys[k].data();
      };

      // The old seed for layer k included every earlier row/s transcript.  Once
      // all COT material is fixed, snapshot one root instead and RO-expand it by
      // the public batch shape and layer index.  Each layer restarts MITCCRH's gid,
      // so the layer separation is load-bearing.
      const block root = fs_coin_pair_group<nP>(
          io, grp, peer, "AGMPC laand batch seed", session_id);
      std::vector<block> layer_seed(B);
      for (int k = 0; k < B; ++k)
        layer_seed[k] =
            RO("emp-ag:laand-layer-seed", session_id)
                .absorb(root)
                .absorb((uint64_t)B)
                .absorb((uint64_t)L)
                .absorb((uint64_t)k)
                .squeeze_block();

      // Keep only the gb-byte T window needed after the bulk s open. phi is a
      // cheap local expression over resident layer material, so each direction
      // recomputes it instead of retaining B*L extra blocks.
      if (bkt_window_.size() < (size_t)BL * gb)
        bkt_window_.resize((size_t)BL * gb);
      std::vector<unsigned char>(&s)[nP + 1] = bkt_s_;
      for (int p = 1; p <= nP; ++p)
        if (s[p].size() < (size_t)BL) s[p].resize((size_t)BL);

#ifdef AG_PROFILE
      uint64_t _tp_hg0 = agmpc_now_ns();
      int64_t _phi0 = io->count();
#endif
      // One task per (lane,direction) walks every layer in canonical order.
      // Do not enqueue per-layer tasks on one NetIO: their writes could interleave.
      const int hg_chunk = 1 << 16;
      int T = std::min(L_lanes, std::max(1, L));
#ifdef EMP_AG_LAAND_LANE_CAP
      T = std::min(T, EMP_AG_LAAND_LANE_CAP);
#endif
      vector<future<void>> rows;
      for (int lane = 0; lane < T; ++lane) {
        auto [lo, hi] = ag::lane_range(lane, L, T);
        NetIOMP<nP> *ch = lane_mesh(lane);
        rows.push_back(pool->enqueue([&, ch, lane, lo, hi]() {
          const int len = hi - lo;
          std::vector<unsigned char> gbuf((size_t)std::min(hg_chunk, len) * gb);
          block pad[16];
          for (int layer = 0; layer < B; ++layer) {
            BlockVec *mac = layer_mac(layer);
            BlockVec *key = layer_key(layer);
            MITCCRH<8> mitc;
            mitc.setS(layer_seed[layer] ^ makeBlock(0, (uint64_t)lane));
            for (int st = lo; st < hi; st += hg_chunk) {
              const int ed = std::min(st + hg_chunk, hi);
              for (int k0 = st; k0 < ed; k0 += 8) {
                const int batch = std::min(8, ed - k0);
                for (int j = 0; j < 8; ++j) {
                  const block kk = (j < batch) ? key[peer][k0 + j] : zero_block;
                  pad[2 * j] = kk;
                  pad[2 * j + 1] = kk ^ Delta;
                }
                mitc.hash_cir<8, 2>(pad);
                for (int j = 0; j < batch; ++j) {
                  const int k = k0 + j;
                  const block phi =
                      (select_mask[getLSB(mac[peer][L + k])] & Delta) ^
                      key[peer][L + k] ^ mac[peer][L + k];
                  const block row = pad[2 * j] ^ pad[2 * j + 1] ^ phi;
                  memcpy(gbuf.data() + (size_t)(k0 + j - st) * gb, &row, gb);
                }
              }
              ch->send_data(peer, gbuf.data(), (size_t)(ed - st) * gb);
            }
          }
          ch->flush(peer);
        }));
        rows.push_back(pool->enqueue([&, ch, lane, lo, hi]() {
          const int len = hi - lo;
          std::vector<unsigned char> gbuf((size_t)std::min(hg_chunk, len) * gb);
          BlockVec wire(std::min(hg_chunk, len));
          block pad[16];
          for (int layer = 0; layer < B; ++layer) {
            BlockVec *mac = layer_mac(layer);
            BlockVec *key = layer_key(layer);
            MITCCRH<8> mitc;
            mitc.setS(layer_seed[layer] ^ makeBlock(0, (uint64_t)lane));
            for (int st = lo; st < hi; st += hg_chunk) {
              const int ed = std::min(st + hg_chunk, hi);
              ch->recv_data(peer, gbuf.data(), (size_t)(ed - st) * gb);
              for (int r = 0; r < ed - st; ++r) {
                wire[r] = zero_block;
                memcpy(&wire[r], gbuf.data() + (size_t)r * gb, gb);
              }
              for (int k0 = st; k0 < ed; k0 += 8) {
                const int batch = std::min(8, ed - k0);
                for (int j = 0; j < 8; ++j) {
                  pad[2 * j] = (j < batch) ? mac[peer][k0 + j] : zero_block;
                  pad[2 * j + 1] =
                      (j < batch) ? key[peer][k0 + j] : zero_block;
                }
                mitc.hash_cir<8, 2>(pad);
                for (int j = 0; j < batch; ++j) {
                  const int k = k0 + j;
                  const block E =
                      pad[2 * j] ^
                      (wire[k - st] & select_mask[getLSB(mac[peer][k])]);
                  block value = pad[2 * j + 1] ^ E ^ key[peer][2 * L + k] ^
                                mac[peer][2 * L + k];
                  const block phi =
                      (select_mask[getLSB(mac[peer][L + k])] & Delta) ^
                      key[peer][L + k] ^ mac[peer][L + k];
                  value ^= phi & select_mask[getLSB(mac[peer][k])];
                  value ^= Delta & select_mask[getLSB(mac[peer][2 * L + k])];
#ifdef EMP_AG_TEST_FLIP_LAAND
                  // Deterministically corrupt ONE late-layer receiver value.
                  // Flipping a shipped row is masked away when its random selector
                  // is zero, which makes an ell==1 adversarial test probabilistic.
                  if (layer == B - 1 && k == 0 && party == 1 && lane == 0)
                    value ^= bit0_mask;
#endif
                  const size_t flat = (size_t)layer * L + k;
                  s[party][flat] = LSB1(value);
                  memcpy(bkt_window_.data() + flat * gb, &value, gb);
                }
              }
            }
          }
        }));
      }
      joinNclean(rows);
#ifdef AG_PROFILE
      g_ag_phi_bytes += (uint64_t)(io->count() - _phi0);
      g_ag_tp_hg_ns += agmpc_now_ns() - _tp_hg0;
      uint64_t _tp_so0 = agmpc_now_ns();
#endif

      // The single dependent flight: open every layer's s vector together.
      lane_fan_out(BL,
        [this, &s](NetIOMP<nP> *ch, int dst, int lo, int hi) {
          ch->send_bool(dst, (const bool *)(s[party].data() + lo), hi - lo);
          ch->flush(dst);
        },
        [&s](NetIOMP<nP> *ch, int src, int lo, int hi) {
          ch->recv_bool(src, (bool *)(s[src].data() + lo), hi - lo);
        });
#ifdef AG_PROFILE
      g_ag_tp_sopen_ns += agmpc_now_ns() - _tp_so0;
      uint64_t _tp_ac0 = agmpc_now_ns();
#endif

      const block dxor = Delta ^ bit0_mask;
      unsigned char delta_bytes[sizeof(block)];
      memcpy(delta_bytes, &Delta, sizeof(block));
      for (int layer = 0; layer < B; ++layer) {
        BlockVec *mac = layer_mac(layer);
        BlockVec *key = layer_key(layer);
        unsigned char *window =
            bkt_window_.data() + (size_t)layer * L * gb;
        ag::parallel_for(pool, L, [&](int lo, int hi) {
          for (int k = lo; k < hi; ++k) {
            const size_t flat = (size_t)layer * L + k;
            const unsigned char opened = s[1][flat] ^ s[2][flat];
            if (opened) {
              if (party == 1)
                mac[peer][2 * L + k] ^= bit0_mask;
              else
                key[1][2 * L + k] ^= dxor;
              for (int j = 0; j < gb; ++j)
                window[(size_t)k * gb + j] ^= delta_bytes[j];
            }
          }
        });
        ag::accumulate_leaky_and_window_bytes(feq, window, L, gb, pool);
      }
#ifdef AG_PROFILE
      g_ag_tp_acheck_ns += agmpc_now_ns() - _tp_ac0;
#endif
    }
  }

  void mint_random_leaky_layer(BlockVec layerMAC[nP + 1],
                               BlockVec layerKEY[nP + 1], int L) {
#ifdef AG_PROFILE
    uint64_t _tp_c0 = agmpc_now_ns();
#endif
    abit.compute(layerMAC, layerKEY, 3 * L);
#ifdef AG_PROFILE
    g_ag_tp_cot_ns += agmpc_now_ns() - _tp_c0;
#endif
  }

  void generate_random_leaky_layer(BlockVec layerMAC[nP + 1],
                                   BlockVec layerKEY[nP + 1], int L, FsHash &feq) {
    mint_random_leaky_layer(layerMAC, layerKEY, L);
    produce_leaky_ands_halfgate(layerMAC, layerKEY, L, feq);
  }

  // The BucketPolicy: the n-party seams behind the shared batched schedule
  // (run_bucketing_core, above). The rep
  // (layer 0, accMAC/accKEY) is filled + leaky'd + stashed by compute_inplace
  // before run_bucketing_core; this policy owns the B-1 per-peer sacrifice buffers,
  // the NetIOMP mesh d-exchange + cross-peer fold, the F_coin (fs_coin_pair /
  // sampleRandom), and the terminal check. nP==2 finalizes with feq_check; nP>=3
  // with the ONE batched RLC zero-sum over the stashed windowed alpha vectors, its
  // challenge phi_cat CONTINUING the post-shift PRG stream (so it rides the same
  // post-commitment coin as the shifts).
  struct BucketPolicy {
    TriplePool *tp;
    BlockVec *accMAC, *accKEY;   // layer 0 (per-peer accumulator slots)
    int L, B, any_peer;
    FsHash *feq;
    std::vector<std::vector<BlockVec>> &sac_macs, &sac_keys;  // tp's reused [B][nP+1]
    std::vector<unsigned char> (&d)[nP + 1];                  // tp's reused (B-1)*L
    FeqSymmetric2PC feq2;                                     // nP==2 only
    std::vector<int> checked_shifts;

    BucketPolicy(TriplePool *tp, BlockVec *accMAC, BlockVec *accKEY, int L, int B,
                 FsHash &feq, int any_peer)
        : tp(tp), accMAC(accMAC), accKEY(accKEY), L(L), B(B), any_peer(any_peer),
          feq(&feq), sac_macs(tp->bkt_sac_macs_), sac_keys(tp->bkt_sac_keys_),
          d(tp->bkt_bd_), checked_shifts(B, 0) {
      if ((int)sac_macs.size() < B) { sac_macs.resize(B); sac_keys.resize(B); }
      for (int p = 1; p <= nP; ++p)
        if (d[p].size() < (size_t)(B - 1) * L) d[p].resize((size_t)(B - 1) * L);
    }

    void ensure_layer_(int k) {
      if ((int)sac_macs[k].size() < nP + 1) {
        sac_macs[k].resize(nP + 1);
        sac_keys[k].resize(nP + 1);
      }
    }

    // Generate all sacrifices. At nP==2, mint every COT first and let the
    // streamlined producer batch layer 0 plus all sacrifices into two flights.
    // The multiparty path keeps its established per-layer FZero/check schedule.
    void generate_layers() {
      if constexpr (nP == 2) {
        for (int k = 1; k < B; ++k) {
          ensure_layer_(k);
          tp->mint_random_leaky_layer(sac_macs[k].data(), sac_keys[k].data(), L);
        }
        tp->produce_leaky_ands_halfgate_batch_2pc(
            accMAC, accKEY, sac_macs, sac_keys, B, L, *feq);
      } else {
        for (int k = 1; k < B; ++k) {
          ensure_layer_(k);
          tp->generate_random_leaky_layer(sac_macs[k].data(), sac_keys[k].data(),
                                          L, *feq);
        }
      }
    }

    // ONE post-commitment coin: FS (nP==2) / sampleRandom (nP>=3).
    block draw_coins() {
      block c;
      if constexpr (nP == 2)
        c = fs_coin_pair_group<nP>(tp->io, tp->grp, any_peer,
                                   "AGMPC bucket", tp->session_id);
      else c = sampleRandom<nP>(tp->io, &tp->prg, tp->pool, tp->party,
                                tp->session_id);
      // Piggyback the half-gate hash seed off this ONE post-commitment coin: a single
      // extra domain-separated squeeze, no new round and no dedicated seed coin. Every
      // party folds the SAME coin here, so the evaluator's half-gate re-hash matches
      // each garbler. Consumed by AGMPCProtocol::draw_and_seed via ctx.mitc.setS.
      tp->halfgate_seed_ = fs_coin("AGMPC halfgate seed", c, zero_block);
      return c;
    }

    // Width-split one layer's folds across the pool, peer-major: each task owns
    // one peer's [i0, i1) output range (streams only that peer's acc/sac arrays),
    // and the d-bit fill rides alongside (it reads the b-region, which no fold
    // writes). One join per layer; layers stay ordered because both RMW the same
    // acc slots. These folds stream hundreds of MB per chunk and were the one
    // serial section left in the bucketing schedule.
    void fold_ranges_(int L_, const std::function<void(int, int, int)> &body) {
      const int nPeers = nP - 1;
      const int per_peer = std::max(1, (int)tp->pool->size() / nPeers);
      const int step = (L_ + per_peer - 1) / per_peer;
      vector<future<void>> fres;
      for (int peer = 1; peer <= nP; ++peer) if (peer != tp->party)
        for (int i0 = 0; i0 < L_; i0 += step) {
          const int i1 = std::min(i0 + step, L_);
          fres.push_back(tp->pool->enqueue([&body, peer, i0, i1]() {
            body(peer, i0, i1);
          }));
        }
      joinNclean(fres);
    }
    void fold_pre_and_dbit(int k, int sh) {
      checked_shifts[k] = sh;
      unsigned char *dk = d[tp->party].data() + (size_t)(k - 1) * L;
      fold_ranges_(L, [this, k, sh, dk](int peer, int i0, int i1) {
        ag::bucket_fold_pre_ac(accMAC[peer].data(), accKEY[peer].data(),
                               sac_macs[k][peer].data(), sac_keys[k][peer].data(), L,
                               sh, i0, i1);
        if (peer == any_peer)
          ag::bucket_fold_dbit(accMAC[any_peer].data(), sac_macs[k][any_peer].data(),
                               L, sh, dk, i0, i1);
      });
    }
    // ONE all-to-all mesh exchange of the full (B-1)*L d-vector + cross-peer fold.
    void exchange_all_d() {
      const int BLd = (B - 1) * L;
#ifdef EMP_AG_LOCAL_TEST_HOOKS
      if (tp->local_open_test_hook)
        tp->local_open_test_hook("bucket-d", d[tp->party].data(), BLd);
#endif
      if constexpr (nP == 2) {
        char Dme[FsHash::DIGEST_SIZE];
        feq->digest(Dme);
        feq2 = feq_prepare_symmetric_2pc(
            &tp->prg, tp->session_id, "emp-ag:laand", tp->party, Dme);
      }
#ifdef AG_PROFILE
      uint64_t _b0 = agmpc_now_ns();
#endif
      // Lane-split the d-vector transfer across the group meshes (§A.3, LANE_CLEAN: a plain
      // bit-vector, no keying, no coin reads it -- the XOR-combine below runs on the ASSEMBLED
      // d[] after the one join, unchanged).
      tp->lane_fan_out(BLd,
        [this](NetIOMP<nP> *ch, int peer, int lo, int hi) {
          ch->send_bool(peer, (const bool *)(d[tp->party].data() + lo), hi - lo);
          if constexpr (nP == 2)
            if (ch == tp->io)
              ch->send_data(peer, feq2.commitment, FsHash::DIGEST_SIZE);
          ch->flush(peer);
        },
        [this](NetIOMP<nP> *ch, int peer, int lo, int hi) {
          ch->recv_bool(peer, (bool *)(d[peer].data() + lo), hi - lo);
          if constexpr (nP == 2)
            if (ch == tp->io)
              ch->recv_data(peer, feq2.peer_commitment, FsHash::DIGEST_SIZE);
        });
      // Local hardening: KRRW's eventual circuit checker is not available
      // to stand-alone GMW/WRK. Authenticate the declared bucket differences
      // before they affect any product share.
      tp->abit.maybe_flush_cot_check();
      const unsigned char *checked_views[nP + 1]{};
      for (int p = 1; p <= nP; ++p) checked_views[p] = d[p].data();
      authenticate_open_claims(tp->io, tp->pool, tp->party, tp->Delta,
          tp->session_id, "triple-bucket-d", tp->checked_open_epoch_++, BLd,
          checked_views,
          [&](int peer, size_t j) {
            const int k = static_cast<int>(j / L) + 1, i = j % L;
            const int src = bucket_fold_src(i, L, checked_shifts[k]);
            return accMAC[peer][L + i] ^ sac_macs[k][peer][L + src];
          },
          [&](int peer, size_t j) {
            const int k = static_cast<int>(j / L) + 1, i = j % L;
            const int src = bucket_fold_src(i, L, checked_shifts[k]);
            return accKEY[peer][L + i] ^ sac_keys[k][peer][L + src];
          });
      // Combined public d = XOR over all parties. Plain range split (NOT
      // fold_ranges_: that is peer-major and would re-apply the whole p-loop
      // once per peer, cancelling the XOR). d[1] ranges are disjoint across
      // tasks; every d[p] is read-only here.
      {
        const int step = (BLd + (int)tp->pool->size() - 1) / (int)tp->pool->size();
        vector<future<void>> cres;
        for (int j0 = 0; j0 < BLd; j0 += step) {
          const int j1 = std::min(j0 + step, BLd);
          cres.push_back(tp->pool->enqueue([this, j0, j1]() {
            for (int p = 2; p <= nP; ++p)
              for (int j = j0; j < j1; ++j)
                d[1][j] = d[1][j] != d[p][j];
          }));
        }
        joinNclean(cres);
      }
#ifdef AG_PROFILE
      g_ag_tp_bkt_ns += agmpc_now_ns() - _b0;
#endif
    }
    void fold_post(int k, int sh) {
      const unsigned char *dk = d[1].data() + (size_t)(k - 1) * L;
      fold_ranges_(L, [this, k, sh, dk](int peer, int i0, int i1) {
        ag::bucket_fold_post(accMAC[peer].data(), accKEY[peer].data(),
                             sac_macs[k][peer].data(), sac_keys[k][peer].data(), L, sh,
                             dk, i0, i1);
      });
    }
    // nP>=3: ONE batched RLC zero-sum over all B layers' stashed windowed alpha
    // vectors (phi_cat continues the post-shift PRG). nP==2: terminal feq_check.
    void finalize_check(PRG &prg) {
      const int party = tp->party;
      if constexpr (nP >= 3) {
        const int BL = (int)tp->alpha_masked_cat_.size();   // = B * L
        BlockVec &phi_cat = tp->bkt_phi_cat_;
        if (phi_cat.size() < (size_t)BL) phi_cat.resize((size_t)BL);
        prg.random_block(phi_cat.data(), BL);
        BlockVec ip[nP + 1];
        for (int i = 0; i <= nP; ++i) ip[i].resize(2);
        // Batched-RLC inner product over all B*L stashed alpha vectors -- the
        // heaviest single-thread compute in the nP>=3 draw+seed critical path.
        // Split into partial (lo,hi) sums and XOR-combine: GF(2^128) addition is
        // XOR and the product sum has no reduction here, so the split is
        // byte-identical to the sequential inner product.
        {
          const int W = std::max(1, (int)tp->pool->size());
          const int width = (BL + W - 1) / W;
          std::vector<block> acc_lo(W, zero_block), acc_hi(W, zero_block);
          vector<future<void>> rres;
          int t = 0;
          for (int a = 0; a < BL; a += width, ++t) {
            const int b = std::min(a + width, BL);
            rres.push_back(tp->pool->enqueue([&, a, b, t]() {
              block part[2];
              vector_inn_prdt_sum_no_red(part, phi_cat.data() + a,
                                         tp->alpha_masked_cat_.data() + a, b - a);
              acc_lo[t] = part[0];
              acc_hi[t] = part[1];
            }));
          }
          joinNclean(rres);
          ip[party][0] = zero_block;
          ip[party][1] = zero_block;
          for (int w = 0; w < t; ++w) {
            ip[party][0] = ip[party][0] ^ acc_lo[w];
            ip[party][1] = ip[party][1] ^ acc_hi[w];
          }
        }
        block uz = zero_block;
        for (block u : tp->alpha_uzero_) uz = uz ^ u;   // all B layers, incl. rep
        ip[party][0] = ip[party][0] ^ uz;
        char dgst[nP + 1][FsHash::DIGEST_SIZE];
        role_bound_commitment(dgst[party], tp->session_id,
                              "emp-ag:laand-alpha", party,
                              ip[party].data(), sizeof(block) * 2);
        vector<future<void>> ares;
        for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
          ares.push_back(tp->pool->enqueue([this, &dgst, peer, party]() {
            tp->io->send_data(peer, dgst[party], FsHash::DIGEST_SIZE);
            tp->io->flush(peer);
          }));
          ares.push_back(tp->pool->enqueue([this, &dgst, peer]() {
            tp->io->recv_data(peer, dgst[peer], FsHash::DIGEST_SIZE);
          }));
        }
        joinNclean(ares);
        char view_digest[FsHash::DIGEST_SIZE];
        FsHash::hash_once(view_digest, dgst[1], nP * FsHash::DIGEST_SIZE);
        char peer_view[nP + 1][FsHash::DIGEST_SIZE];
        vector<future<bool>> ares2;
        for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
          ares.push_back(tp->pool->enqueue([this, &ip, &view_digest, peer, party]() {
            tp->io->send_data(peer, ip[party].data(), sizeof(block) * 2);
            tp->io->send_data(peer, view_digest, FsHash::DIGEST_SIZE);
            tp->io->flush(peer);
          }));
          ares2.push_back(tp->pool->enqueue(
              [this, &ip, &dgst, &peer_view, &view_digest, peer]() -> bool {
            tp->io->recv_data(peer, ip[peer].data(), sizeof(block) * 2);
            tp->io->recv_data(peer, peer_view[peer], FsHash::DIGEST_SIZE);
            char chk[FsHash::DIGEST_SIZE];
            role_bound_commitment(chk, tp->session_id,
                                  "emp-ag:laand-alpha", peer,
                                  ip[peer].data(), sizeof(block) * 2);
            return memcmp(chk, dgst[peer], FsHash::DIGEST_SIZE) != 0 ||
                   memcmp(view_digest, peer_view[peer],
                          FsHash::DIGEST_SIZE) != 0;
          }));
        }
        joinNclean(ares);
        expecting(!joinNcleanCheat(ares2), "LaAND alpha: commit-open mismatch");
        for (int i = 2; i <= nP; ++i)
          xorBlocks_arr(ip[1].data(), ip[1].data(), ip[i].data(), 2);
        expecting(cmpBlock(&ip[1][0], &zero_block, 1) &&
                  cmpBlock(&ip[1][1], &zero_block, 1),
                  "LaAND alpha: Sigma alpha != 0");
      }
      if constexpr (nP == 2) {
        feq_finish_symmetric_2pc(
            tp->send_ch(any_peer), tp->recv_ch(any_peer), tp->session_id,
            "emp-ag:laand", party, feq2,
            "AGMPC LaAND F_eq: leaky-AND check failed");
      }
    }
  };

  // Size the per-peer accumulators for a length-L Beaver batch and run the int-
  // overflow guard, BEFORE the draw walk: DrawBeaverPass (passes/draw_beaver.h) captures the
  // Beaver operands directly into acc's a/b regions ([0, L) and [L, 2L)) as it
  // visits each AND, so the arrays must be at their 3L size when the walk starts.
  void prepare_beaver_acc(int length) {
    if (length == 0) return;
    const int B = get_bucket_size(length);
    // All bucketing arithmetic scales L by up to max(3,B) in int (3*L sacrifice
    // mints, (B-1)*L d-vectors, B*L RLC stash). Reject a flush too large for that to
    // fit int rather than silently wrapping -- callers split work into chunks long
    // before this bound for memory reasons anyway.
    expecting((int64_t)std::max(B, 3) * length <= (int64_t)INT_MAX,
              "triple_pool: flush too large (bucketed size overflows int) -- split into chunks");
    for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
      if (bkt_accMAC_[peer].size() < (size_t)3 * length) {
        bkt_accMAC_[peer].resize((size_t)3 * length);
        bkt_accKEY_[peer].resize((size_t)3 * length);
      }
    }
  }

  // KRRW-style function-dependent AND share generation. The representative layer
  // is the real pair of input masks, already sitting in acc's a/b regions (the
  // draw walk captured them there; prepare_beaver_acc sized the arrays), with only
  // r freshly minted; B-1 sacrifice layers are generated and bucketed one at a
  // time. a_bit/b_bit are the original per-AND share bits (bit0 of the slot-0
  // MAC), needed by the final Beaver x/y open after the folds clobber acc.
  void compute_inplace(const unsigned char *a_bit, const unsigned char *b_bit,
                       int length, AShareBundleVec<nP> &out_sigma) {
    out_sigma.clear();
    if (length == 0) return;
    int B = get_bucket_size(length);
    int L = length;
    int any_peer = (party == 1) ? 2 : 1;

#ifdef AG_PROFILE
    g_ag_tp_cot_ns = g_ag_tp_hg_ns = g_ag_tp_sopen_ns = 0;
    g_ag_tp_acheck_ns = g_ag_tp_bkt_ns = g_ag_tp_open_ns = 0;
#endif

    BlockVec(&accMAC)[nP + 1] = bkt_accMAC_;
    BlockVec(&accKEY)[nP + 1] = bkt_accKEY_;

    // Retain the authenticated originals for the final checked Beaver open.
    // Bucketing changes acc.a/acc.b; share bits alone cannot authenticate
    // the original-minus-folded differences.
    for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
      originalMAC_[peer].assign(accMAC[peer].begin(), accMAC[peer].begin() + 2 * L);
      originalKEY_[peer].assign(accKEY[peer].begin(), accKEY[peer].begin() + 2 * L);
    }

    BlockVec(&rMAC)[nP + 1] = bkt_rMAC_;
    BlockVec(&rKEY)[nP + 1] = bkt_rKEY_;
#ifdef AG_PROFILE
    uint64_t _tp_c0 = agmpc_now_ns();
#endif
    abit.compute(rMAC, rKEY, L);
#ifdef AG_PROFILE
    g_ag_tp_cot_ns += agmpc_now_ns() - _tp_c0;
#endif
    for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
      memcpy(accMAC[peer].data() + 2 * L, rMAC[peer].data(), L * sizeof(block));
      memcpy(accKEY[peer].data() + 2 * L, rKEY[peer].data(), L * sizeof(block));
    }
    FsHash feq;  // running F_eq transcript over all layers' gb-windows (nP==2 only)
    alpha_masked_cat_.clear();   // nP>=3 batched-RLC stash, fresh per bucket
    alpha_uzero_.clear();
    if constexpr (nP >= 3)
      produce_leaky_ands_halfgate(accMAC, accKEY, L, feq); // rep layer 0 (stashes)

    // Bucketing via the shared batched schedule (run_bucketing_core, above):
    // nP==2 batches all B layer rows/s; nP>=3 generates the B-1 sacrifices,
    // then ONE post-commitment coin drives all shifts, followed by one mesh d exchange,
    // then the terminal check (nP==2 feq_check / nP>=3 batched RLC over the stashed
    // windowed alpha vectors). The rep (layer 0, produced + stashed just above) and
    // the B-1 sacrifices all feed the check. See BucketPolicy for the n-party seams.
    BucketPolicy pol(this, accMAC, accKEY, L, B, feq, any_peer);
    ag::run_bucketing_core(pol, B, L);

#ifdef AG_PROFILE
    uint64_t _tp_o0 = agmpc_now_ns();
#endif
    std::vector<unsigned char>(&x)[nP + 1] = bkt_x_;
    std::vector<unsigned char>(&y)[nP + 1] = bkt_y_;
    for (int p = 1; p <= nP; ++p) {
      if (x[p].size() < (size_t)L) x[p].resize((size_t)L);
      if (y[p].size() < (size_t)L) y[p].resize((size_t)L);
    }
    ag::parallel_for(pool, L, [&](int lo, int hi) {
      for (int i = lo; i < hi; ++i) {
        x[party][i] = (unsigned char)(a_bit[i] ^ getLSB(accMAC[any_peer][i]));
        y[party][i] = (unsigned char)(b_bit[i] ^ getLSB(accMAC[any_peer][L + i]));
      }
    });
#ifdef EMP_AG_LOCAL_TEST_HOOKS
    if (local_open_test_hook) local_open_test_hook("final-x", x[party].data(), L);
#endif
    // Lane-split the Beaver x/y open (§A.3, LANE_CLEAN; AFTER run_bucketing_core so no
    // coin/check reads it -- the combine below runs on the assembled x[]/y[]). Both x and y
    // ride the same slice on each mesh, so P1's window cannot desync.
    lane_fan_out(L,
      [this, &x, &y](NetIOMP<nP> *ch, int peer, int lo, int hi) {
        ch->send_bool(peer, (const bool *)(x[party].data() + lo), hi - lo);
        ch->send_bool(peer, (const bool *)(y[party].data() + lo), hi - lo);
        ch->flush(peer);
      },
      [&x, &y](NetIOMP<nP> *ch, int peer, int lo, int hi) {
        ch->recv_bool(peer, (bool *)(x[peer].data() + lo), hi - lo);
        ch->recv_bool(peer, (bool *)(y[peer].data() + lo), hi - lo);
      });
    // Authenticate BOTH Beaver opening vectors before accepting their fold.
    // Domain/epoch and the canonical all-party view are checked as well.
    abit.maybe_flush_cot_check();
    std::vector<std::vector<unsigned char>> checked_xy(nP + 1);
    const unsigned char *checked_views[nP + 1]{};
    for (int p = 1; p <= nP; ++p) {
      checked_xy[p].resize(2 * L);
      std::copy_n(x[p].data(), L, checked_xy[p].data());
      std::copy_n(y[p].data(), L, checked_xy[p].data() + L);
      checked_views[p] = checked_xy[p].data();
    }
    authenticate_open_claims(io, pool, party, Delta, session_id,
        "triple-final-xy", checked_open_epoch_++, size_t(2) * L, checked_views,
        [&](int peer, size_t j) { return originalMAC_[peer][j] ^ accMAC[peer][j]; },
        [&](int peer, size_t j) { return originalKEY_[peer][j] ^ accKEY[peer][j]; });
    // x[1]/y[1] combine (disjoint i ranges; each d[p] read-only).
    ag::parallel_for(pool, L, [&](int lo, int hi) {
      for (int p = 2; p <= nP; ++p)
        for (int i = lo; i < hi; ++i) {
          x[1][i] = x[1][i] ^ x[p][i];
          y[1][i] = y[1][i] ^ y[p][i];
        }
    });

    out_sigma.resize(L);
    block dxor = Delta ^ bit0_mask;
    // Branchless Beaver fold, width-split: out_sigma[i] is disjoint per i, and
    // every read (accMAC/accKEY/x/y) is done. This is the heaviest of the
    // final-open loops (full AShareBundle arithmetic per element).
    ag::parallel_for(pool, L, [&](int lo, int hi) {
      for (int i = lo; i < hi; ++i) {
        bool xb = x[1][i], yb = y[1][i];
        // mask each addend by select_mask instead of branching on the two random
        // public selectors (byte-identical to the if(xb)/if(yb) form).
        const block mx = select_mask[xb], my = select_mask[yb];
        AShareBundle<nP> &sb = out_sigma[i];
        for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
          int slot = peer_slot(party, peer);
          sb.mac(slot) = accMAC[peer][2 * L + i] ^ (accMAC[peer][L + i] & mx) ^ (accMAC[peer][i] & my);
          sb.key(slot) = accKEY[peer][2 * L + i] ^ (accKEY[peer][L + i] & mx) ^ (accKEY[peer][i] & my);
        }
        if (xb && yb) {
          if (party != 1) {
            sb.key(peer_slot(party, 1)) = sb.key(peer_slot(party, 1)) ^ dxor;
          } else {
            for (int slot = 0; slot < AShareBundle<nP>::K; ++slot)
              sb.mac(slot) = sb.mac(slot) ^ bit0_mask;
          }
        }
      }
    });
#ifdef AG_PROFILE
    g_ag_tp_open_ns += agmpc_now_ns() - _tp_o0;
    if (party == 1) {
      auto ms = [](uint64_t ns) { return (double)ns / 1.0e6; };
      uint64_t sum = g_ag_tp_cot_ns + g_ag_tp_hg_ns + g_ag_tp_sopen_ns +
                     g_ag_tp_acheck_ns + g_ag_tp_bkt_ns + g_ag_tp_open_ns;
      std::printf("[ag-tp]   -- inplace_triples breakdown (B=%d layers, L=%d) --\n", B, L);
      std::printf("[ag-tp]     COT/aShare draw  %9.1f ms\n", ms(g_ag_tp_cot_ns));
      std::printf("[ag-tp]     half-gate phi    %9.1f ms\n", ms(g_ag_tp_hg_ns));
      std::printf("[ag-tp]     s-open           %9.1f ms\n", ms(g_ag_tp_sopen_ns));
      std::printf("[ag-tp]     alpha check      %9.1f ms\n", ms(g_ag_tp_acheck_ns));
      std::printf("[ag-tp]     bucket layers    %9.1f ms\n", ms(g_ag_tp_bkt_ns));
      std::printf("[ag-tp]     final open       %9.1f ms\n", ms(g_ag_tp_open_ns));
      std::printf("[ag-tp]     sub-sum          %9.1f ms\n", ms(sum));
    }
#endif
  }
};

}  // namespace emp::ag
#endif // EMP_AG_BACKEND_PREPROC_TRIPLE_POOL_H__
