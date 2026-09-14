// src/protocol/sign.h — Pi_MLDSA Sign, single-shot (T = 1), SLOT-RESTRICTED
// boundary test (threshold-mldsa-slot.pdf, Fig. 1).
//
// Window recap (unchanged from the two-round design). With H^z = beta - c*s,
// H^w = beta + c*e (so 0 <= H <= 2*beta) and U^z = gamma1 - y, U^w = gamma2 - w0,
//   accept_i  <=>  2*beta < U_i + H_i < S_i,   S = 2*gamma1 resp. 2*gamma2,
// and the pre-challenge producer emits M_i = R_i + shift(U_i, S_i) mod 2^RING_K
// such that T_i = M_i + p_i = shift + H has its top bit set exactly on accept.
//
// SLOT IDEA: a coefficient can only fail when U_i is within 2*beta of an
// endpoint (d_i = 1); interior coefficients accept for EVERY reachable H. The
// expected number of edge coefficients is ~2.9 of m >= 2048, so the online test
// runs on NS = 64 slots holding only the edge coefficients:
//   Phase 1 (pre-challenge, Boolean)  producers -> (M_i, d_i); a running one-hot
//     of cnt_i = #edges before i gives u_{j,i} = d_i AND [cnt_i = j] at one AND
//     per (i, j); an order-preserving compaction network carries M into the
//     slots (M^_j), empty slots get L = 2^(RING_K-1) so they accept; ovf =
//     [more than NS edges] is opened before the challenge.
//   Phase 2 (fixed, garbled before c)  K_j = M^_j + d^_j, r = AND_j K_j[top],
//     Z_out = r * bin(rho).
// The response no longer rides on the boundary sum: rho is a fresh one-time pad
// and V = (gamma1+beta) - y - c*s + rho is opened in flight 1 (safe, rho is
// uniform and used once); Z_out releases rho only when r = 1, so a rejected
// attempt reveals nothing but the bit.
//
// Online: one checked BDOZ ring opening carrying (delta_H, V); the routing
// d^_j[l] = XOR_{i: delta_H,i[l]=1} u_{j,i} is F2-linear with PUBLIC coefficients
// (free XOR on retained wires); then the prepared circuit and one reveal.
#ifndef MLDSA_SIGN_H
#define MLDSA_SIGN_H

#include "circuit.h"
#include "dealer.h"
#include "keygen.h"
#include "prepsign.h"
#include "ref.h"
#include <emp-ag/emp-ag.h>
#include <emp-tool/circuits/frontend/frontend.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace mldsa {

constexpr int NS = 64;                                    // slots (paper Table 1)
constexpr int M_TOTAL = Y_COEFF_COUNT + COEFF_COUNT;      // response + low-bit positions
constexpr int ceil_log2(int v) {
  int r = 0;
  while ((1 << r) < v)
    ++r;
  return r;
}
constexpr int LOGM = ceil_log2(M_TOTAL);
constexpr int DW = LOGM;                                  // shift amount width
constexpr int HALF = 1 << (RING_K - 1);
constexpr int64_t SPAN_Z = 2 * (int64_t)GAMMA1;           // U^z in [0, SPAN_Z)
constexpr int64_t SPAN_W = 2 * (int64_t)G2;               // U^w in [0, SPAN_W)
constexpr int64_t C_LOW = 2 * BETA + 1;
constexpr int64_t C_J = HALF - 2 * BETA - 1;
constexpr int64_t Z_OFF = (int64_t)GAMMA1 + BETA;         // V - rho = Z_OFF - z (mod q)
static_assert(SPAN_Z == (int64_t{1} << Y_WIDTH), "U^z must be exactly Y_WIDTH bits");
static_assert(SPAN_W < (int64_t{1} << OW0), "U^w must fit OW0 bits");

constexpr int MHW = NS * RING_K;                  // M^ slots
constexpr int DHW = NS * RING_K;                  // routed residues
constexpr int RHOW = Y_COEFF_COUNT * L;           // bin23(rho), the F_q response pad
constexpr int OUTW = 1 + RHOW;                    // r || Z_out

template <class Ctx> using MHV = emp::BitVec_T<Ctx, MHW>;
template <class Ctx> using DHV = emp::BitVec_T<Ctx, DHW>;
template <class Ctx> using RHOV = emp::BitVec_T<Ctx, RHOW>;
template <class Ctx> using OUTV = emp::BitVec_T<Ctx, OUTW>;

// (c, z, h) with explicit bot markers: z_ok=false => (c,bot,bot) [rejected or
// slot overflow]; h_ok=false => (c,z,bot). All public — this IS the signature.
struct Signature {
  std::vector<int32_t> c; // challenge polynomial, N coefficients in {-1,0,1}
  bool z_ok = false;
  std::vector<int32_t> z; // N*ELL coefficients, |z| < GAMMA1-BETA (iff z_ok)
  bool h_ok = false;
  std::vector<uint8_t> h; // N*K hint bits (iff h_ok)
};

// acc[N] += c * x over R_q on pairwise-authenticated field shares. Local.
template <int nP>
inline void cpoly_mul_fq(const std::vector<std::pair<int, int>>& c_nz, const FqShare<nP>* x,
                         FqShare<nP>* acc) {
  for (const auto& [deg, sgn0] : c_nz)
    for (int i = 0; i < N; ++i) {
      int d = i - deg, sgn = sgn0;
      if (d < 0) {
        d += N;
        sgn = -sgn;
      }
      acc[i] = acc[i] + (sgn > 0 ? x[d] : x[d] * (uint32_t)(Q - 1));
    }
}

// ---- local ring arithmetic --------------------------------------------------

// acc[N] += c * x over Z_{2^RING_W}[X]/(X^N+1) on ring shares; c sparse in
// {-1,0,1} given as (degree, sign) pairs. Purely local (linear).
template <int nP>
inline void cpoly_mul_ring(const std::vector<std::pair<int, int>>& c_nz, const RingShare<nP>* x,
                           RingShare<nP>* acc) {
  for (const auto& [deg, sgn0] : c_nz)
    for (int i = 0; i < N; ++i) {
      int d = i - deg, sgn = sgn0;
      if (d < 0) { // negacyclic wrap: x^N = -1
        d += N;
        sgn = -sgn;
      }
      acc[i] = sgn > 0 ? acc[i] + x[d] : acc[i] - x[d];
    }
}

// ---- Phase 1 (pre-challenge Boolean) ----------------------------------------

// From U (UW bits, unsigned, < span) and the eDaBit bits R (RING_K bits):
//   lower = [U <= 2beta], upper = [U >= span - 2beta], edge = lower ^ upper,
//   J = HALF + edge*(U - HALF) + C_J*lower - (span mod 2^RING_K)*upper,
//   M = R + J mod 2^RING_K.
// Costs 2 compares + RING_K + (RING_K-1) AND; `edge` is free for the caller.
template <int UW, class Ctx>
inline void window_producer(Ctx& ctx, const Bit_T<Ctx>* U, const Bit_T<Ctx>* R, int64_t span,
                            Bit_T<Ctx>& edge_out, Bit_T<Ctx>* M) {
  using Bit = Bit_T<Ctx>;
  const Bit lower = !ge_const_u<UW>(ctx, U, (uint64_t)C_LOW);
  const Bit upper = ge_const_u<UW>(ctx, U, (uint64_t)(span - 2 * BETA));
  const Bit edge = lower ^ upper; // the two cases are exclusive (2beta < span/2)
  edge_out = edge;

  Bit t[RING_K];
  for (int k = 0; k < RING_K; ++k)
    t[k] = edge & (k == RING_K - 1 ? !U[k] : U[k]);
  t[RING_K - 1] = !t[RING_K - 1]; // + HALF
  const int64_t NSPAN = (-(span & (int64_t)RING_KMASK)) & (int64_t)RING_KMASK;
  Bit J[RING_K], carry = Bit::constant(ctx, false);
  for (int k = 0; k < RING_K; ++k) {
    const bool cj = ((C_J >> k) & 1) != 0, ns = ((NSPAN >> k) & 1) != 0;
    if (cj || ns)
      fa(J[k], carry, t[k], cj && ns ? edge : (cj ? lower : upper), carry);
    else
      ha(J[k], carry, t[k], carry);
  }
  add_ripple<RING_K>(ctx, M, R, J, false);
}

// Order-preserving compaction. Item i carries payload A[i] (VW bits, zero when
// empty) and shift S[i] (DW bits, zero when empty); it must end at i - S[i].
// LSB first: at stage b every item whose remaining shift has bit b set moves
// left by 2^b, realised as a conditional swap of (p, p+2^b) driven by the
// routed shift bit. Because the destination map is monotone, the cell a real
// item moves into is always free at that moment, and sweeping p upward lets a
// whole chain shift in one stage. Cost: one AND per swapped bit,
// sum_b (m - 2^b) * (VW + DW - b).
template <int VW, class Ctx>
inline void compact(std::vector<std::array<Bit_T<Ctx>, VW>>& A,
                    std::vector<std::array<Bit_T<Ctx>, DW>>& S, int m) {
  using Bit = Bit_T<Ctx>;
  for (int b = 0; b < LOGM; ++b) {
    const int step = 1 << b;
    for (int p = 0; p + step < m; ++p) {
      const Bit ctrl = S[(size_t)p + step][(size_t)b];
      for (int k = 0; k < VW; ++k) {
        const Bit t = ctrl & (A[(size_t)p][(size_t)k] ^ A[(size_t)p + step][(size_t)k]);
        A[(size_t)p][(size_t)k] = A[(size_t)p][(size_t)k] ^ t;
        A[(size_t)p + step][(size_t)k] = A[(size_t)p + step][(size_t)k] ^ t;
      }
      for (int k = b; k < DW; ++k) { // bits below b are already consumed
        const Bit t = ctrl & (S[(size_t)p][(size_t)k] ^ S[(size_t)p + step][(size_t)k]);
        S[(size_t)p][(size_t)k] = S[(size_t)p][(size_t)k] ^ t;
        S[(size_t)p + step][(size_t)k] = S[(size_t)p + step][(size_t)k] ^ t;
      }
    }
  }
}

// ---- Phase 2: the fixed post-challenge circuit -------------------------------

template <class Ctx>
inline OUTV<Ctx> phase2_body(Ctx& ctx, const MHV<Ctx>& Mh, const DHV<Ctx>& Dh,
                             const RHOV<Ctx>& Rho) {
  using Bit = Bit_T<Ctx>;
  std::vector<Bit> ok((size_t)NS);
  Bit a[RING_K], b[RING_K], K[RING_K];
  for (int j = 0; j < NS; ++j) {
    for (int k = 0; k < RING_K; ++k) {
      a[k] = Mh[j * RING_K + k];
      b[k] = Dh[j * RING_K + k];
    }
    add_ripple<RING_K>(ctx, K, a, b, false);
    ok[(size_t)j] = K[RING_K - 1];
  }
  const Bit r = and_tree(ok.data(), NS);
  std::vector<Bit> outb((size_t)OUTW);
  outb[0] = r;
  for (int t = 0; t < RHOW; ++t)
    outb[(size_t)(1 + t)] = r & Rho[t]; // Z_out = r * bin(rho)
  return OUTV<Ctx>::from_bit_values(ctx, outb.data());
}

// Compiled once per process: the circuit depends on the parameter set only.
inline const auto& phase2_circuit() {
  using emp::RecordCtx;
  static const auto c =
      emp::frontend::compile<MHV<RecordCtx>, DHV<RecordCtx>, RHOV<RecordCtx>>(
          [](RecordCtx& ctx, MHV<RecordCtx> Mh, DHV<RecordCtx> Dh, RHOV<RecordCtx> Rho) {
            return phase2_body(ctx, Mh, Dh, Rho);
          });
  return c;
}

// No-op default for sign()'s after_prepsign hook (see below).
struct NoopHook {
  void operator()() const {}
};

// ---- Sign -------------------------------------------------------------------
// `after_prepsign` fires at the offline/online boundary (after F_PrepSign, the
// eDaBits, Phase 1 and the prepared garbling). `online_rounds_out` receives the
// number of online synchronization barriers; `g2_ands_out` the AND count of the
// prepared post-challenge circuit.
template <int nP, class Hook = NoopHook>
inline Signature sign(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer,
                      KeyPair<nP>& kp, const std::vector<uint32_t>& msg,
                      Hook after_prepsign = Hook{}, int64_t* online_rounds_out = nullptr,
                      int64_t* g2_ands_out = nullptr) {
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using Bit = emp::Bit_T<Ctx>;
  using Wire = typename Ctx::Wire;
  using UR = emp::UInt_T<Ctx, RING_K>;

  Signature sig;

#ifdef TEST
  auto lap_start = emp::clock_start();
  auto timer_lap = [&](const char* label) {
    if (party == 1)
      std::printf("  [timer] %-20s %8.1f ms\n", label, emp::time_from(lap_start) / 1000.0);
    lap_start = emp::clock_start();
  };
#endif

  // ---- offline ---------------------------------------------------------------
  // 1a. F_PrepSign -> (<w0>_2 wires, <y> shares, w1 public), T = 1.
  std::vector<emp::UInt_T<Ctx, OW0>> w0_2;
  SharePair<nP, Y_WIDTH, false> y;
  std::vector<uint32_t> w1;
  const uint64_t and_a0 = sess.num_and();
  prepsign<nP>(sess, party, dealer, kp.A, w0_2, y, w1);
  const uint64_t and_a1 = sess.num_and();
#ifdef TEST
  timer_lap("prepsign");
#endif

  // (R_H, bin_nu) per position over the ring, and the F_q response pad
  // (rho, bin23(rho)) per response coefficient -- paper Fig. 1, step 1.
  RingEdabits<nP> rh = ring_edabits<nP, RING_K>(sess, party, dealer, M_TOTAL);
  SharePair<nP, L, true> rho = Fq_edabits<nP, L>(sess, party, dealer, Y_COEFF_COUNT);

  // Phase 1: producers, running one-hot of cnt, shift prefix, compaction.
  std::vector<std::array<Bit, RING_K>> Marr((size_t)M_TOTAL);
  std::vector<std::array<Bit, DW>> Sarr((size_t)M_TOTAL);
  std::vector<Bit> edge((size_t)M_TOTAL);
  {
    Bit U[Y_WIDTH > OW0 ? Y_WIDTH : OW0], R[RING_K];
    for (int idx = 0; idx < M_TOTAL; ++idx) {
      const auto r = sess.template adopt_authenticated_input<UR>(&rh.two_share[(size_t)RING_K * idx]);
      for (int k = 0; k < RING_K; ++k)
        R[k] = r[k];
      if (idx < Y_COEFF_COUNT) {
        // <y>_2 stores v = y - 1 (two's complement); U = gamma1 - y = ~v with
        // the top bit restored -- pure rewiring.
        const auto v = sess.template adopt_authenticated_input<emp::UInt_T<Ctx, Y_WIDTH>>(
            &y.two_share[(size_t)Y_WIDTH * idx]);
        for (int k = 0; k < Y_WIDTH; ++k)
          U[k] = k < Y_WIDTH - 1 ? !v[k] : v[k];
        window_producer<Y_WIDTH>(sess.ctx(), U, R, SPAN_Z, edge[(size_t)idx],
                                 Marr[(size_t)idx].data());
      } else {
        Bit w0b[OW0], nw0[OW0];
        unpack_bits<OW0>(w0_2[(size_t)(idx - Y_COEFF_COUNT)], w0b);
        for (int k = 0; k < OW0; ++k)
          nw0[k] = !w0b[k];
        add_const_ripple<OW0>(sess.ctx(), U, nw0, (uint64_t)G2, true); // U = gamma2 - w0
        window_producer<OW0>(sess.ctx(), U, R, SPAN_W, edge[(size_t)idx],
                             Marr[(size_t)idx].data());
      }
    }
  }
#ifdef TEST
  timer_lap("phase1 producers");
#endif

  // u[j][i] = edge_i AND [cnt_i = j], carried by a running one-hot h of cnt.
  // h_{i+1}[j] = h_i[j] ^ u[j-1][i] ^ u[j][i] is free once u is formed, so the
  // whole decoder costs one AND per (i, j) -- m*(NS+1) in total.
  std::vector<std::vector<Bit>> u((size_t)NS, std::vector<Bit>((size_t)M_TOTAL));
  {
    std::vector<Bit> h((size_t)NS + 1, Bit::constant(sess.ctx(), false));
    h[0] = Bit::constant(sess.ctx(), true);
    std::vector<Bit> ucol((size_t)NS + 1);
    Bit shift[DW];
    for (int k = 0; k < DW; ++k)
      shift[k] = Bit::constant(sess.ctx(), false);
    for (int i = 0; i < M_TOTAL; ++i) {
      for (int j = 0; j <= NS; ++j)
        ucol[(size_t)j] = edge[(size_t)i] & h[(size_t)j];
      for (int j = 0; j < NS; ++j)
        u[(size_t)j][(size_t)i] = ucol[(size_t)j];
      for (int j = NS; j >= 0; --j)
        h[(size_t)j] = h[(size_t)j] ^ (j > 0 ? ucol[(size_t)j - 1] : Bit::constant(sess.ctx(), false)) ^
                       ucol[(size_t)j];
      // route only edge items: shift = edge ? (#non-edges so far) : 0
      for (int k = 0; k < DW; ++k)
        Sarr[(size_t)i][(size_t)k] = edge[(size_t)i] & shift[k];
      add_bit_ripple<DW>(sess.ctx(), shift, shift, !edge[(size_t)i]);
    }
    // h is one-hot while cnt <= NS and all-zero once it overflows (the 1 is
    // shifted out of the NS slot), so ovf is a free XOR reduction.
    Bit acc = h[0];
    for (int j = 1; j <= NS; ++j)
      acc = acc ^ h[(size_t)j];
    const auto ovf = sess.reveal(!acc, emp::PUBLIC);
    emp::expecting(ovf.has_value(), "Sign: ovf reveal failed");
    if (ovf.value()) {
      std::fprintf(stderr, "Sign: slot overflow (>%d edge coefficients), retire attempt\n", NS);
      return sig; // (c, bot, bot) -- happens with probability ~2^-207
    }
  }
#ifdef TEST
  timer_lap("phase1 decoder");
#endif

  compact<RING_K>(Marr, Sarr, M_TOTAL);
#ifdef TEST
  timer_lap("phase1 compaction");
#endif

  // Slot j = compacted position j. Positions beyond the edge count still hold
  // an interior coefficient's producer output, so an empty slot is MUXed to
  // L = 2^(RING_K-1), which makes K_j[top] = 1 (accept). NS*RING_K AND.
  std::vector<Wire> mhw((size_t)MHW);
  for (int j = 0; j < NS; ++j) {
    Bit empty = Bit::constant(sess.ctx(), true);
    for (int i = 0; i < M_TOTAL; ++i)
      empty = empty ^ u[(size_t)j][(size_t)i];
    for (int k = 0; k < RING_K; ++k) {
      const Bit m = Marr[(size_t)j][(size_t)k];
      const Bit bit = k == RING_K - 1 ? (m ^ (empty & !m)) : (m & !empty);
      bit.pack_wires(&mhw[(size_t)j * RING_K + k]);
    }
  }
  const MHV<Ctx> Mh_v = MHV<Ctx>::from_wires(sess.ctx(), mhw.data());

  // Boolean half of the response pad, adopted once.
  std::vector<Wire> rhow((size_t)RHOW);
  {
    using URHO = emp::UInt_T<Ctx, L>;
    for (int i = 0; i < Y_COEFF_COUNT; ++i)
      sess.template adopt_authenticated_input<URHO>(&rho.two_share[(size_t)L * i])
          .pack_wires(&rhow[(size_t)i * L]);
  }
  const RHOV<Ctx> Rho_v = RHOV<Ctx>::from_wires(sess.ctx(), rhow.data());
  sess.checkpoint(); // materialise Phase 1: u, M^ and rho become soldering sources
  if (party == 1) // G1 split: prepsign vs producers + decoder + compaction
    std::printf("G1SPLIT prepsign=%llu phase1=%llu\n", (unsigned long long)(and_a1 - and_a0),
                (unsigned long long)(sess.num_and() - and_a1));

  const uint64_t ands_before_prepare = sess.num_and();
  auto prepared = sess.prepare(phase2_circuit());
  if (g2_ands_out)
    *g2_ands_out = (int64_t)(sess.num_and() - ands_before_prepare);
#ifdef TEST
  timer_lap("prepare phase 2");
#endif
  after_prepsign(); // offline/online boundary -- everything below needs msg

  // ---- online ----------------------------------------------------------------
  // 1b. mu = H(tr, m), c = H(mu, w1).
  std::vector<uint32_t> mu_in(kp.tr.begin(), kp.tr.end());
  mu_in.insert(mu_in.end(), msg.begin(), msg.end());
  std::array<uint32_t, 8> mu;
  h256_stub(mu_in, mu);
  std::vector<uint32_t> c_in(mu.begin(), mu.end());
  c_in.insert(c_in.end(), w1.begin(), w1.end());
  std::array<uint32_t, 8> c_seed;
  h256_stub(c_in, c_seed);
  sample_in_ball_stub(c_seed, sig.c);

  std::vector<std::pair<int, int>> c_nz; // (degree, sign), TAU entries
  for (int j = 0; j < N; ++j)
    if (sig.c[j])
      c_nz.push_back({j, sig.c[j]});

  // Flight 1 (paper step 5): one opening carrying both domains --
  //   [[delta_H]]_D = (beta*1 - c[[s]]_D, beta*1 + c[[e]]_D) - [[R_H]]_D   (ring)
  //   [[V]]_q       = (gamma1+beta)*1 - [[y]]_q - c[[s]]_q + [[rho]]_q     (field)
  // Both halves carry pairwise (BDOZ) MACs verified locally by the receiver, so
  // the whole opening is a single exchange -- the paper's atomic checked Open.
  std::vector<RingShare<nP>> dh((size_t)M_TOTAL);
  std::vector<FqShare<nP>> vbar_sh((size_t)Y_COEFF_COUNT);
  {
    std::vector<RingShare<nP>> cs((size_t)Y_COEFF_COUNT), ce((size_t)COEFF_COUNT);
    for (int p = 0; p < ELL; ++p)
      cpoly_mul_ring(c_nz, &kp.s.ring_share[(size_t)p * N], &cs[(size_t)p * N]);
    for (int p = 0; p < K; ++p)
      cpoly_mul_ring(c_nz, &kp.e.ring_share[(size_t)p * N], &ce[(size_t)p * N]);
    for (int i = 0; i < Y_COEFF_COUNT; ++i)
      dh[(size_t)i] = add_public_ring(RingShare<nP>{} - cs[(size_t)i], (uint64_t)BETA, party,
                                      dealer.my_alpha_r) -
                      rh.ring_share[(size_t)i];
    for (int i = 0; i < COEFF_COUNT; ++i)
      dh[(size_t)Y_COEFF_COUNT + i] =
          add_public_ring(ce[(size_t)i], (uint64_t)BETA, party, dealer.my_alpha_r) -
          rh.ring_share[(size_t)Y_COEFF_COUNT + i];

    std::vector<FqShare<nP>> csq((size_t)Y_COEFF_COUNT);
    for (int p = 0; p < ELL; ++p)
      cpoly_mul_fq(c_nz, &kp.s.fq_share[(size_t)p * N], &csq[(size_t)p * N]);
    for (int i = 0; i < Y_COEFF_COUNT; ++i) {
      FqShare<nP> v = y.fq_share[(size_t)i] * (uint32_t)(Q - 1);
      v = v + csq[(size_t)i] * (uint32_t)(Q - 1) + rho.fq_share[(size_t)i];
      vbar_sh[(size_t)i] = add_public_fq(v, (uint32_t)(Z_OFF % Q), party, dealer.my_alpha_f);
    }
  }
#ifdef TAMPER_C
  if (party == 1)
    dh[0].val.w[0] ^= 1;
#endif
  std::vector<RingVal> opened;
  std::vector<uint32_t> vbar;
  open_ring_and_fq(sess.io(), party, dealer.my_alpha_r, dealer.my_alpha_f, dh, vbar_sh, opened,
                   vbar);

  // Routing: d^_j[l] = XOR over the public index set { i : delta_H,i[l] = 1 }.
  // Free XOR on retained wires; no PublicBits, no per-coefficient online work.
  std::vector<Wire> dhw((size_t)DHW);
  for (int l = 0; l < RING_K; ++l) {
    std::vector<int> idx;
    idx.reserve((size_t)M_TOTAL);
    for (int i = 0; i < M_TOTAL; ++i)
      if ((opened[(size_t)i].w[0] >> l) & 1)
        idx.push_back(i);
    for (int j = 0; j < NS; ++j) {
      Bit acc = Bit::constant(sess.ctx(), false);
      for (int i : idx)
        acc = acc ^ u[(size_t)j][(size_t)i];
      acc.pack_wires(&dhw[(size_t)j * RING_K + l]);
    }
  }
  const DHV<Ctx> Dh_v = DHV<Ctx>::from_wires(sess.ctx(), dhw.data());
  sess.checkpoint(); // settle the XOR-only chunk so d^ can be soldered

  const OUTV<Ctx> out = sess.run(std::move(prepared), Mh_v, Dh_v, Rho_v);
  const auto released = sess.reveal(out, 1); // Open((r, Z_out), P1)
  if (online_rounds_out)
    *online_rounds_out = 4; // flight 1, routing checkpoint, run, reveal
  if (party != 1)
    return sig;
  emp::expecting(released.has_value(), "Sign: online reveal failed");
#ifdef TEST
  timer_lap("flight1 + routing + run");
#endif

  const auto& ob = released.value();
  const bool accept = ob[0];

#ifdef TEST
  // Plaintext oracle: recompute r0 / z from the opened e, s (keygen) and w, y
  // (prepsign); the circuit's accept bit and (if accepted) z must agree.
  const auto centered = [](uint32_t x) {
    return x > (uint32_t)Q / 2 ? (int32_t)x - Q : (int32_t)x;
  };
  bool expect_accept = true;
  for (int p = 0; p < K; ++p) {
    int64_t e_row[N], ce_row[N] = {0};
    for (int i = 0; i < N; ++i)
      e_row[i] = centered(kp.opened_e[p * N + i]);
    cpoly_mul_acc(sig.c, e_row, ce_row);
    for (int i = 0; i < N; ++i) {
      int64_t d = ref_lowbits((int32_t)test_opened_w[p * N + i], G2) - ce_row[i];
      if (d > G2)
        d -= 2 * (int64_t)G2;
      else if (d <= -G2)
        d += 2 * (int64_t)G2;
      if (d >= G2 - BETA || d <= -(G2 - BETA))
        expect_accept = false;
    }
  }
  std::vector<int64_t> z_plain(Y_COEFF_COUNT);
  for (int p = 0; p < ELL; ++p) {
    int64_t s_row[N], cs_row[N] = {0};
    for (int i = 0; i < N; ++i)
      s_row[i] = centered(kp.opened_s[p * N + i]);
    cpoly_mul_acc(sig.c, s_row, cs_row);
    for (int i = 0; i < N; ++i) {
      const int idx = p * N + i;
      z_plain[idx] = centered(test_opened_y[idx]) + cs_row[i];
      if (z_plain[idx] >= GAMMA1 - BETA || z_plain[idx] <= -(GAMMA1 - BETA))
        expect_accept = false;
    }
  }
  if (party == 1)
    std::printf("  Sign: circuit accept=%d, plaintext accept=%d\n", (int)accept,
                (int)expect_accept);
  emp::expecting(accept == expect_accept, "Sign test: accept bit disagrees with plaintext");
  if (!accept)
    for (int t = 1; t < OUTW; ++t)
      emp::expecting(!ob[(size_t)t], "Sign test: rejected run released pad material");
#endif

  if (!accept)
    return sig; // (c, bot, bot)

  // z = (gamma1+beta)*1 - ((V - Z_out) mod q)   (paper step 7).
  sig.z.resize(Y_COEFF_COUNT);
  for (int i = 0; i < Y_COEFF_COUNT; ++i) {
    uint32_t zout = 0;
    for (int k = 0; k < L; ++k)
      zout |= (uint32_t)ob[(size_t)(1 + i * L + k)] << k;
    const uint32_t v = fq_sub(vbar[(size_t)i], zout);
    sig.z[i] = (int32_t)(Z_OFF - (int64_t)v);
#ifdef TEST
    emp::expecting(sig.z[i] == (int32_t)z_plain[i], "Sign test: opened z != plaintext y + c*s");
#endif
  }
  sig.z_ok = true;

  // 4.+5. All-public tail: h = MakeHint(-c t0, Az - ct + c t0, 2g2), then
  // ||c t0||_inf >= g2 or HW(h) > omega -> (c, z, bot).
  std::vector<uint32_t> z_modq(Y_COEFF_COUNT);
  for (int i = 0; i < Y_COEFF_COUNT; ++i)
    z_modq[i] = sig.z[i] >= 0 ? (uint32_t)sig.z[i] : (uint32_t)(sig.z[i] + Q);
  const std::vector<uint32_t> Az = matvec_pub_negacyclic(kp.A, z_modq, K, ELL);

  sig.h.assign(COEFF_COUNT, 0);
  int hint_weight = 0;
  int64_t ct0_max = 0;
  for (int p = 0; p < K; ++p) {
    int64_t t_row[N], t0_row[N], ct[N] = {0}, ct0[N] = {0};
    for (int i = 0; i < N; ++i) {
      t_row[i] = (int64_t)kp.t[p * N + i];
      t0_row[i] = (int64_t)kp.t0[p * N + i];
    }
    cpoly_mul_acc(sig.c, t_row, ct);
    cpoly_mul_acc(sig.c, t0_row, ct0);
    for (int i = 0; i < N; ++i) {
      const int idx = p * N + i;
      const int64_t az = (int64_t)Az[idx];
      const uint32_t uu = modq_i64(az - ct[i] + ct0[i]); // Az - ct + c*t0
      const uint8_t hi = (uint8_t)ref_makehint(-ct0[i], (int32_t)uu, G2);
      sig.h[idx] = hi;
      hint_weight += hi;
      const int64_t mag = ct0[i] < 0 ? -ct0[i] : ct0[i];
      ct0_max = mag > ct0_max ? mag : ct0_max;
    }
  }
  sig.h_ok = ct0_max < G2 && hint_weight <= OMEGA;
#ifdef TEST
  if (party == 1)
    std::printf("  Sign: ||c*t0||_inf=%lld (g2=%d), HW(h)=%d (omega=%d)\n",
                (long long)ct0_max, G2, hint_weight, OMEGA);
#endif
  if (!sig.h_ok) {
    sig.h.clear();
    return sig; // (c, z, bot)
  }

#ifdef TEST
  // End-to-end oracle: run the ML-DSA verifier. w1' = UseHint(h, Az - c*2^d*t1)
  // must reproduce w1 and therefore c. The paper's variant carries e_w slack
  // (Az - ct = w - e_w - c*e), so report rather than abort on a mismatch.
  std::vector<uint32_t> w1p(COEFF_COUNT);
  for (int p = 0; p < K; ++p) {
    int64_t t1s_row[N], ct1s[N] = {0};
    for (int i = 0; i < N; ++i)
      t1s_row[i] = (int64_t)kp.t1[p * N + i] << POW2ROUND_D;
    cpoly_mul_acc(sig.c, t1s_row, ct1s);
    for (int i = 0; i < N; ++i) {
      const int idx = p * N + i;
      const uint32_t wapp = modq_i64((int64_t)Az[idx] - ct1s[i]);
      w1p[idx] = (uint32_t)ref_usehint(sig.h[idx], (int32_t)wapp, G2);
    }
  }
  std::vector<uint32_t> v_in(mu.begin(), mu.end());
  v_in.insert(v_in.end(), w1p.begin(), w1p.end());
  std::array<uint32_t, 8> v_seed;
  h256_stub(v_in, v_seed);
  std::vector<int32_t> c_verify;
  sample_in_ball_stub(v_seed, c_verify);
  int w1_diff = 0;
  for (int i = 0; i < COEFF_COUNT; ++i)
    w1_diff += w1p[i] != w1[i];
  if (party == 1)
    std::printf("  Sign verify: w1 mismatches %d / %d, challenge %s\n", w1_diff, COEFF_COUNT,
                c_verify == sig.c ? "MATCHES" : "DIFFERS");
  timer_lap("hint + verify");
#endif

  return sig;
}

} // namespace mldsa
#endif // MLDSA_SIGN_H
