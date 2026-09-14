// src/infra/bdoz.h — the arithmetic layer: pairwise-authenticated (BDOZ) shares over F_q
// and the widened power-of-two ring, with receiver-verified single-flight openings.
//
// The world in FRONT of the GC. Every authenticated value -- over F_q and over
// the widened power-of-two ring -- carries PAIRWISE (BDOZ) MACs: for each
// ordered pair (j, i) party P_j holds a one-time key on P_i's share and P_i
// holds M = K + alpha_j * share. An opening is therefore verified locally by
// the receiver and costs a single flight, in both domains. F_q MACs are
// repeated FQ_MAC_REP times because one 23-bit MAC only gives 2^-23.
// FakeDealer (bottom of this file) FAKES the offline from a shared seed --
// replace it (dealer.h) with a real offline; keep the
// share layout and the checked opens.
#ifndef MLDSA_BDOZ_H
#define MLDSA_BDOZ_H

#include "ref.h"                     // mldsa::Q, fq_* field ops
#include <emp-ag/backend/netmp.h>    // emp::NetIOMP (vendored WRK/GMW mesh)
#include <emp-tool/emp-tool.h>       // emp::expecting, emp::Hash, emp::block

#include <array>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

namespace mldsa {

using emp::expecting;
using emp::NetIOMP;

// ---- authenticated shares over the ring Z_{2^(nu+sigma)} ---------------------
// Widened power-of-two ring for the post-challenge window (paper Sec. 3.1:
// BDOZ-style pairwise authentication with the SPDZ2k widened representation).
// The semantic value lives in Z_{2^RING_K} (nu = b+1 bits); the extra
// RING_SIGMA bits carry the authentication slack: an alteration of a low-nu
// residue has 2-adic valuation at most nu-1, so a forgery succeeds with
// probability about 2^-RING_SIGMA. MLDSA_SIGMA is a build parameter
// (make SIGMA=...); the papers instantiate sigma = 128.
constexpr int RING_K = 2 * BETA < 256 ? 9 : 10; // nu = b + 1, with L = 2^b > 2*BETA
static_assert(2 * BETA < (1 << (RING_K - 1)), "RING_K window too narrow for 2*beta");
#ifndef MLDSA_SIGMA
#define MLDSA_SIGMA 128
#endif
constexpr int RING_SIGMA = MLDSA_SIGMA;
constexpr int RING_W = RING_K + RING_SIGMA; // nu + sigma
constexpr uint64_t RING_KMASK = (uint64_t{1} << RING_K) - 1;
constexpr int RING_WORDS = (RING_W + 63) / 64;
static_assert(RING_SIGMA >= 8 && RING_WORDS <= 3, "sigma must keep nu+sigma within 192 bits");

// Fixed-width unsigned integer mod 2^RING_W, little-endian 64-bit words.
struct RingVal {
  uint64_t w[RING_WORDS] = {};
};
inline void ring_reduce(RingVal& a) {
  constexpr int top = RING_W - 64 * (RING_WORDS - 1);
  if constexpr (top < 64)
    a.w[RING_WORDS - 1] &= (uint64_t{1} << top) - 1;
}
inline RingVal ring_add(const RingVal& a, const RingVal& b) {
  RingVal r;
  unsigned __int128 c = 0;
  for (int i = 0; i < RING_WORDS; ++i) {
    c += (unsigned __int128)a.w[i] + b.w[i];
    r.w[i] = (uint64_t)c;
    c >>= 64;
  }
  ring_reduce(r);
  return r;
}
inline RingVal ring_sub(const RingVal& a, const RingVal& b) {
  RingVal r;
  unsigned __int128 borrow = 0;
  for (int i = 0; i < RING_WORDS; ++i) {
    const unsigned __int128 d = (unsigned __int128)a.w[i] - b.w[i] - borrow;
    r.w[i] = (uint64_t)d;
    borrow = (d >> 64) ? 1 : 0;
  }
  ring_reduce(r);
  return r;
}
// Schoolbook low-half product: only the low RING_W bits are defined.
inline RingVal ring_mul(const RingVal& a, const RingVal& b) {
  uint64_t acc[RING_WORDS + 1] = {};
  for (int i = 0; i < RING_WORDS; ++i) {
    uint64_t carry = 0;
    for (int j = 0; i + j < RING_WORDS; ++j) {
      const unsigned __int128 t =
          (unsigned __int128)a.w[i] * b.w[j] + acc[i + j] + carry;
      acc[i + j] = (uint64_t)t;
      carry = (uint64_t)(t >> 64);
    }
    if (i + RING_WORDS <= RING_WORDS)
      acc[RING_WORDS] += carry;
  }
  RingVal r;
  for (int i = 0; i < RING_WORDS; ++i)
    r.w[i] = acc[i];
  ring_reduce(r);
  return r;
}
inline RingVal ring_from_u64(uint64_t v) {
  RingVal r;
  r.w[0] = v;
  ring_reduce(r);
  return r;
}
inline RingVal ring_from_i64(int64_t v) { // two's complement lift
  RingVal r = ring_from_u64((uint64_t)(v < 0 ? -v : v));
  return v < 0 ? ring_sub(RingVal{}, r) : r;
}
inline uint64_t ring_low(const RingVal& a) { return a.w[0]; } // low 64 bits

// ---- wire format: RING_W bits per value, bit-packed (paper: m(nu+sigma) bits)
inline size_t ring_packed_bytes(size_t count) { return (count * (size_t)RING_W + 7) / 8; }
inline std::vector<uint8_t> ring_pack(const std::vector<RingVal>& v) {
  std::vector<uint8_t> out(ring_packed_bytes(v.size()), 0);
  size_t bit = 0;
  for (const RingVal& x : v)
    for (int k = 0; k < RING_W; ++k, ++bit)
      if ((x.w[k >> 6] >> (k & 63)) & 1)
        out[bit >> 3] |= (uint8_t)(1u << (bit & 7));
  return out;
}
inline std::vector<RingVal> ring_unpack(const std::vector<uint8_t>& in, size_t count) {
  std::vector<RingVal> out(count);
  size_t bit = 0;
  for (size_t i = 0; i < count; ++i)
    for (int k = 0; k < RING_W; ++k, ++bit)
      if ((in[bit >> 3] >> (bit & 7)) & 1)
        out[i].w[k >> 6] |= uint64_t{1} << (k & 63);
  return out;
}

// Pairwise (BDOZ) authentication: for every ordered pair (j, i), P_j holds a
// one-time key K_{j,i} on P_i's share and P_i holds M_{j,i} = K_{j,i} + alpha_j * x_i,
// with alpha_j P_j's private slope reused across values. Everything is linear,
// and an opening is verified by the RECEIVER alone -- one flight, no second open.
template <int nP> struct RingShare {
  RingVal val;                       // my additive share of x mod 2^RING_W
  std::array<RingVal, nP + 1> mac{}; // mac[j]: MAC on my share under P_j's key
  std::array<RingVal, nP + 1> key{}; // key[i]: my key on P_i's share
};
template <int nP> inline RingShare<nP> operator+(const RingShare<nP>& a, const RingShare<nP>& b) {
  RingShare<nP> r;
  r.val = ring_add(a.val, b.val);
  for (int p = 1; p <= nP; ++p) {
    r.mac[p] = ring_add(a.mac[p], b.mac[p]);
    r.key[p] = ring_add(a.key[p], b.key[p]);
  }
  return r;
}
template <int nP> inline RingShare<nP> operator-(const RingShare<nP>& a, const RingShare<nP>& b) {
  RingShare<nP> r;
  r.val = ring_sub(a.val, b.val);
  for (int p = 1; p <= nP; ++p) {
    r.mac[p] = ring_sub(a.mac[p], b.mac[p]);
    r.key[p] = ring_sub(a.key[p], b.key[p]);
  }
  return r;
}
template <int nP> inline RingShare<nP> operator*(const RingShare<nP>& x, uint64_t c) {
  const RingVal cv = ring_from_u64(c);
  RingShare<nP> r;
  r.val = ring_mul(x.val, cv);
  for (int p = 1; p <= nP; ++p) {
    r.mac[p] = ring_mul(x.mac[p], cv);
    r.key[p] = ring_mul(x.key[p], cv);
  }
  return r;
}
// x + k for public k: party 1 adds k to its share; every other party lowers its
// key on party 1's share by alpha_me * k so M = K + alpha * x_1 still holds.
template <int nP>
inline RingShare<nP> add_public_ring(RingShare<nP> x, uint64_t k, int party,
                                     const RingVal& alpha_me) {
  if (party == 1)
    x.val = ring_add(x.val, ring_from_u64(k));
  else
    x.key[1] = ring_sub(x.key[1], ring_mul(alpha_me, ring_from_u64(k)));
  return x;
}

// ---- pairwise (BDOZ) authentication over F_q --------------------------------
// A single 23-bit MAC gives only 2^-23 soundness, so the slopes are repeated
// FQ_MAC_REP times with independent keys (paper Sec. 3.1: "the realization must
// amplify authentication ... through suitable extension-field or repeated
// authentication"). Verification is local at the receiver, so an opening needs
// no second flight -- that is what lets V ride in flight 1.
constexpr int FQ_MAC_REP = (MLDSA_SIGMA + L - 1) / L; // R with 2^-(23R) soundness
using FqMac = std::array<uint32_t, FQ_MAC_REP>;
inline FqMac fq_mac_add(const FqMac& a, const FqMac& b) {
  FqMac r{};
  for (int i = 0; i < FQ_MAC_REP; ++i)
    r[i] = fq_add(a[i], b[i]);
  return r;
}
inline FqMac fq_mac_mul(const FqMac& a, uint32_t c) {
  FqMac r{};
  for (int i = 0; i < FQ_MAC_REP; ++i)
    r[i] = fq_mul(a[i], c);
  return r;
}

template <int nP> struct FqShare {
  uint32_t val = 0;                  // my additive share of x mod q
  std::array<FqMac, nP + 1> mac{};   // mac[j] = key_{j,me} + alpha_j * val_me
  std::array<FqMac, nP + 1> key{};   // key[i]: my keys on P_i's share
};
template <int nP> inline FqShare<nP> operator+(const FqShare<nP>& a, const FqShare<nP>& b) {
  FqShare<nP> r;
  r.val = fq_add(a.val, b.val);
  for (int p = 1; p <= nP; ++p) {
    r.mac[p] = fq_mac_add(a.mac[p], b.mac[p]);
    r.key[p] = fq_mac_add(a.key[p], b.key[p]);
  }
  return r;
}
template <int nP> inline FqShare<nP> operator*(const FqShare<nP>& x, uint32_t c) {
  FqShare<nP> r;
  r.val = fq_mul(x.val, c);
  for (int p = 1; p <= nP; ++p) {
    r.mac[p] = fq_mac_mul(x.mac[p], c);
    r.key[p] = fq_mac_mul(x.key[p], c);
  }
  return r;
}
template <int nP>
inline FqShare<nP> add_public_fq(FqShare<nP> x, uint32_t k, int party, const FqMac& alpha_me) {
  if (party == 1)
    x.val = fq_add(x.val, k);
  else
    for (int i = 0; i < FQ_MAC_REP; ++i)
      x.key[1][i] = fq_sub(x.key[1][i], fq_mul(alpha_me[i], k));
  return x;
}

// Checked pairwise opening over F_q, one flight: shares plus one digest of the
// MAC vector per ordered pair, verified locally by the receiver.
template <int nP>
inline std::vector<uint32_t> open_fq_checked(NetIOMP<nP>& io, int party, const FqMac& alpha_me,
                                             const std::vector<FqShare<nP>>& x) {
  const size_t len = x.size();
  std::vector<uint32_t> mine(len);
  for (size_t t = 0; t < len; ++t)
    mine[t] = x[t].val;
  std::vector<uint32_t> macs(len * FQ_MAC_REP);
  auto digest_of = [&](int from, int to, char* out) {
    emp::Hash h;
    const uint32_t ids[2] = {(uint32_t)from, (uint32_t)to};
    h.put(ids, sizeof(ids));
    h.put(macs.data(), (int64_t)macs.size() * sizeof(uint32_t));
    h.digest(out);
  };
  char dg[emp::Hash::DIGEST_SIZE];
  for (int j = 1; j <= nP; ++j)
    if (j != party) {
      for (size_t t = 0; t < len; ++t)
        for (int r = 0; r < FQ_MAC_REP; ++r)
          macs[t * FQ_MAC_REP + r] = x[t].mac[j][r];
      digest_of(party, j, dg);
      io.send_data(j, mine.data(), len * sizeof(uint32_t));
      io.send_data(j, dg, sizeof(dg));
    }
  io.flush();
  std::vector<uint32_t> acc = mine, buf(len);
  char peer_dg[emp::Hash::DIGEST_SIZE];
  for (int i = 1; i <= nP; ++i)
    if (i != party) {
      io.recv_data(i, buf.data(), len * sizeof(uint32_t));
      io.recv_data(i, peer_dg, sizeof(peer_dg));
      for (size_t t = 0; t < len; ++t)
        for (int r = 0; r < FQ_MAC_REP; ++r)
          macs[t * FQ_MAC_REP + r] = fq_add(x[t].key[i][r], fq_mul(alpha_me[r], buf[t]));
      digest_of(i, party, dg);
      expecting(std::memcmp(dg, peer_dg, sizeof(dg)) == 0,
                "BDOZ F_q opening failed: peer share does not match its MAC");
      for (size_t t = 0; t < len; ++t)
        acc[t] = fq_add(acc[t], buf[t]);
    }
  return acc;
}

// ---- public matrix x shared vector over R_q --------------------------------
// acc += A * x in R_q^rows with R_q = F_Q[X]/(X^N+1); A is a public rows x cols
// matrix of polynomials (flattened, N coefficients each), x a shared cols-vector.
// Call with acc preloaded with the additive term (e.g. <e>_q), length rows*N.
// Purely local: FqShare is linear under public-constant multiplication.
template <int nP>
inline void matvec_negacyclic(const std::vector<uint32_t>& A, const std::vector<FqShare<nP>>& x,
                              std::vector<FqShare<nP>>& acc, int rows, int cols) {
  expecting((int)A.size() == rows * cols * N, "matvec_negacyclic: A size");
  expecting((int)x.size() == cols * N, "matvec_negacyclic: x size");
  expecting((int)acc.size() == rows * N, "matvec_negacyclic: acc size");
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      for (int a_deg = 0; a_deg < N; ++a_deg) {
        const uint32_t a = A[((size_t)row * cols + col) * N + a_deg];
        for (int x_deg = 0; x_deg < N; ++x_deg) {
          const int degree = a_deg + x_deg;
          const int out_deg = degree < N ? degree : degree - N;
          // X^N = -1: coefficients that wrap flip sign.
          const uint32_t a_mod_q = degree < N || a == 0 ? a : (uint32_t)Q - a;
          acc[row * N + out_deg] = acc[row * N + out_deg] + x[col * N + x_deg] * a_mod_q;
        }
      }
    }
  }
}

// ---- transport opens -------------------------------------------------------
// One all-to-all round: send my whole vector to every peer, then combine in
// each peer's. Runs BEFORE any GC op and is symmetric + fully drained, so it is
// safe on the session's GC mesh.

// Checked BDOZ ring opening in ONE flight: bit-packed widened shares plus one
// digest of the MAC vector per ordered pair, verified locally by the receiver.
template <int nP>
inline std::vector<RingVal> open_ring_checked(NetIOMP<nP>& io, int party, const RingVal& alpha_me,
                                              const std::vector<RingShare<nP>>& x) {
  const size_t len = x.size();
  std::vector<RingVal> mine(len), tmp(len);
  for (size_t t = 0; t < len; ++t)
    mine[t] = x[t].val;
  const std::vector<uint8_t> mine_packed = ring_pack(mine);
  auto digest_of = [&](int from, int to, const std::vector<RingVal>& macs, char* out) {
    const std::vector<uint8_t> packed = ring_pack(macs);
    emp::Hash h;
    const uint32_t ids[2] = {(uint32_t)from, (uint32_t)to};
    h.put(ids, sizeof(ids));
    h.put(packed.data(), (int64_t)packed.size());
    h.digest(out);
  };
  char dg[emp::Hash::DIGEST_SIZE];
  for (int j = 1; j <= nP; ++j)
    if (j != party) {
      for (size_t t = 0; t < len; ++t)
        tmp[t] = x[t].mac[j];
      digest_of(party, j, tmp, dg);
      io.send_data(j, mine_packed.data(), mine_packed.size());
      io.send_data(j, dg, sizeof(dg));
    }
  io.flush();
  std::vector<RingVal> acc = mine;
  std::vector<uint8_t> buf(mine_packed.size());
  char peer_dg[emp::Hash::DIGEST_SIZE];
  for (int i = 1; i <= nP; ++i)
    if (i != party) {
      io.recv_data(i, buf.data(), buf.size());
      io.recv_data(i, peer_dg, sizeof(peer_dg));
      const std::vector<RingVal> peer = ring_unpack(buf, len);
      for (size_t t = 0; t < len; ++t)
        tmp[t] = ring_add(x[t].key[i], ring_mul(alpha_me, peer[t]));
      digest_of(i, party, tmp, dg);
      expecting(std::memcmp(dg, peer_dg, sizeof(dg)) == 0,
                "BDOZ ring opening failed: peer share does not match its MAC");
      for (size_t t = 0; t < len; ++t)
        acc[t] = ring_add(acc[t], peer[t]);
    }
  return acc;
}

// One-flight opening carrying BOTH domains, as the slot protocol's flight 1
// does. Each half is authenticated pairwise (BDOZ) and verified locally by the
// receiver from its own keys, so the whole thing is a single exchange.
template <int nP>
inline void open_ring_and_fq(NetIOMP<nP>& io, int party, const RingVal& alpha_me,
                             const FqMac& fq_alpha_me, const std::vector<RingShare<nP>>& rx,
                             const std::vector<FqShare<nP>>& fx, std::vector<RingVal>& ring_out,
                             std::vector<uint32_t>& fq_out) {
  const size_t len = rx.size(), flen = fx.size();
  std::vector<uint32_t> fq(flen);
  for (size_t t = 0; t < flen; ++t)
    fq[t] = fx[t].val;
  std::vector<RingVal> mine(len), tmp(len);
  for (size_t t = 0; t < len; ++t)
    mine[t] = rx[t].val;
  const std::vector<uint8_t> mine_packed = ring_pack(mine);
  std::vector<uint32_t> fmac(flen * FQ_MAC_REP);
  auto digest_of = [&](int from, int to, const std::vector<RingVal>& macs, char* out) {
    const std::vector<uint8_t> packed = ring_pack(macs);
    emp::Hash h;
    const uint32_t ids[2] = {(uint32_t)from, (uint32_t)to};
    h.put(ids, sizeof(ids));
    h.put(packed.data(), (int64_t)packed.size());
    h.put(fmac.data(), (int64_t)fmac.size() * sizeof(uint32_t));
    h.digest(out);
  };
  char dg[emp::Hash::DIGEST_SIZE];
  for (int j = 1; j <= nP; ++j)
    if (j != party) {
      for (size_t t = 0; t < len; ++t)
        tmp[t] = rx[t].mac[j];
      for (size_t t = 0; t < flen; ++t)
        for (int i = 0; i < FQ_MAC_REP; ++i)
          fmac[t * FQ_MAC_REP + i] = fx[t].mac[j][i];
      digest_of(party, j, tmp, dg);
      io.send_data(j, mine_packed.data(), mine_packed.size());
      io.send_data(j, fq.data(), flen * sizeof(uint32_t));
      io.send_data(j, dg, sizeof(dg));
    }
  io.flush();
  ring_out = mine;
  fq_out = fq;
  std::vector<uint8_t> buf(mine_packed.size());
  std::vector<uint32_t> fbuf(flen);
  char peer_dg[emp::Hash::DIGEST_SIZE];
  for (int i = 1; i <= nP; ++i)
    if (i != party) {
      io.recv_data(i, buf.data(), buf.size());
      io.recv_data(i, fbuf.data(), flen * sizeof(uint32_t));
      io.recv_data(i, peer_dg, sizeof(peer_dg));
      const std::vector<RingVal> peer = ring_unpack(buf, len);
      for (size_t t = 0; t < len; ++t)
        tmp[t] = ring_add(rx[t].key[i], ring_mul(alpha_me, peer[t]));
      for (size_t t = 0; t < flen; ++t)
        for (int r = 0; r < FQ_MAC_REP; ++r)
          fmac[t * FQ_MAC_REP + r] = fq_add(fx[t].key[i][r], fq_mul(fq_alpha_me[r], fbuf[t]));
      digest_of(i, party, tmp, dg);
      expecting(std::memcmp(dg, peer_dg, sizeof(dg)) == 0,
                "BDOZ opening failed: peer share does not match its MAC");
      for (size_t t = 0; t < len; ++t)
        ring_out[t] = ring_add(ring_out[t], peer[t]);
      for (size_t t = 0; t < flen; ++t)
        fq_out[t] = fq_add(fq_out[t], fbuf[t]);
    }
}

// XOR-open a GF(2) sharing: reconstruct  ⊕_p share_p.  Used to open the drawn r
// (demo cheat) so a consistent ⟨r⟩_q can be built.
template <int nP>
inline std::vector<uint32_t> open_xor(NetIOMP<nP>& io, int party,
                                      const std::vector<uint32_t>& mine) {
  const int len = (int)mine.size();
  std::vector<uint32_t> acc = mine;
  for (int p = 1; p <= nP; ++p)
    if (p != party)
      io.send_data(p, mine.data(), (size_t)len * sizeof(uint32_t));
  io.flush();
  std::vector<uint32_t> buf(len);
  for (int p = 1; p <= nP; ++p)
    if (p != party) {
      io.recv_data(p, buf.data(), (size_t)len * sizeof(uint32_t));
      for (int i = 0; i < len; ++i)
        acc[i] ^= buf[i];
    }
  return acc;
}

} // namespace mldsa
#endif // MLDSA_BDOZ_H
