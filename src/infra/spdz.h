// src/infra/spdz.h — the arithmetic "front" of Pi_PrepSign (demo SPDZ / edaBit layer).
//
// The world in FRONT of the GC: authenticated additive shares over F_Q
// (AuthShare = a value share of x plus a MAC share of α·x under a global key α),
// the transport opens that reconstruct them, and the SPDZ MACCheck that gates a
// public c.  The dealer stubs (arith_split / auth_split) FAKE F_edabits + the
// SPDZ dealer from a shared seed — replace them (and the r-open cheat in main)
// with your real F_edabits / SPDZ engine; keep the opens and the MACCheck.
//
// Templated on party count nP: deduced from NetIOMP<nP> for the opens / MACCheck;
// pass explicitly to the dealer stubs, e.g. arith_split<nP>(rng, v).
#ifndef MLDSA_SPDZ_H
#define MLDSA_SPDZ_H

#include "ref.h"           // mldsa::Q, fq_* field ops
#include <emp-ag/emp-ag.h> // emp::NetIOMP, emp::expecting

#include <array>
#include <cstdint>
#include <random>
#include <vector>

namespace mldsa {

using emp::expecting;
using emp::NetIOMP;

// SPDZ authenticated additive share of x over F_Q: value share of x + MAC share
// of α·x.  Linear ops (add, add-public) are local on both fields.
struct AuthShare {
  uint32_t val = 0; // additive part of x
  uint32_t mac = 0; // additive part of α·x
};
inline AuthShare operator+(const AuthShare& a, const AuthShare& b) {
  return {fq_add(a.val, b.val), fq_add(a.mac, b.mac)}; // linear: both fields add locally
}
inline AuthShare operator*(const AuthShare& x, uint32_t c) {
  return {fq_mul(x.val, c), fq_mul(x.mac, c)};
}

// ---- public matrix x shared vector over R_q --------------------------------
// acc += A * x in R_q^rows with R_q = F_Q[X]/(X^N+1); A is a public rows x cols
// matrix of polynomials (flattened, N coefficients each), x a shared cols-vector.
// Call with acc preloaded with the additive term (e.g. <e>_q), length rows*N.
// Purely local: AuthShare is linear under public-constant multiplication.
inline void matvec_negacyclic(const std::vector<uint32_t>& A, const std::vector<AuthShare>& x,
                              std::vector<AuthShare>& acc, int rows, int cols) {
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

// Reconstruct  Σ_p share_p  (mod Q).
template <int nP>
inline std::vector<uint32_t> open_additive_modq(NetIOMP<nP>& io, int party,
                                                const std::vector<uint32_t>& mine) {
  int len = mine.size();
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
        acc[i] = fq_add(acc[i], buf[i]);
    }
  return acc;
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

// SPDZ batched MACCheck. After opening public c[i], verify  Σ_p mac_i^p == α·c_i
// for all i. Per party: σ_i = mac_i − α^me·c_i; a public random combination
// s = Σ_i χ_i·σ_i is opened across parties and must be 0, because
//   Σ_p s^p = Σ_i χ_i ( Σ_p mac_i^p − (Σ_p α^p)·c_i ) = Σ_i χ_i(α·c_i − α·c_i) = 0.
// A wrong opened c (or tampered mac) makes it nonzero -> abort. THIS is where the
// authenticated shares of w/r "do their work". (DEMO: χ from a fixed public coin
// and σ opened directly; real SPDZ draws χ AFTER commitments and commit-opens σ —
// rush-secure. Same mechanism.)
template <int nP>
inline void spdz_maccheck(NetIOMP<nP>& io, int party, uint32_t my_alpha,
                          const std::vector<uint32_t>& my_mac, const std::vector<uint32_t>& c) {
  expecting(my_mac.size() == c.size(), "spdz_maccheck: my_mac and c same length");
  int len = my_mac.size();
  std::mt19937 coin(0xC0FFEEu); // public coin, same at every party (DEMO)
  uint32_t s = 0;
  for (int i = 0; i < len; ++i) {
    const uint32_t chi = coin() % Q;
    const uint32_t sigma = fq_sub(my_mac[i], fq_mul(my_alpha, c[i])); // mac_i − α·c_i
    s = fq_add(s, fq_mul(chi, sigma));
  }
  const uint32_t sum = open_additive_modq(io, party, std::vector<uint32_t>{s})[0];
  expecting(sum == 0, "SPDZ MACCheck failed: opened c is not authenticated (alpha*c != sum mac)");
}

// ==== DEMO edaBit / SPDZ dealer (replace with F_edabits / real SPDZ) =========
// Deterministic splits from a shared seed: every party COMPUTES the whole split
// but only ever USES its own slot [party] — the visibility a real ⟨·⟩ share gives.
// NOT how real secret sharing works; a self-consistent, checkable stand-in.
template <int nP> inline std::array<uint32_t, nP + 1> arith_split(std::mt19937& rng, uint32_t v) {
  std::array<uint32_t, nP + 1> s{};
  uint32_t acc = 0;
  for (int p = 1; p < nP; ++p) {
    s[p] = rng() % Q;
    acc = fq_add(acc, s[p]);
  }
  s[nP] = fq_sub(v, acc); // last share absorbs the rest
  return s;
}
// SPDZ authenticated split: value shares sum to v, MAC shares sum to α·v.
template <int nP>
inline std::array<AuthShare, nP + 1> auth_split(std::mt19937& rng, uint32_t v, uint32_t alpha) {
  const auto vs = arith_split<nP>(rng, v);                // Σ = v
  const auto ms = arith_split<nP>(rng, fq_mul(alpha, v)); // Σ = α·v
  std::array<AuthShare, nP + 1> s{};
  for (int p = 1; p <= nP; ++p)
    s[p] = {vs[p], ms[p]};
  return s;
}

// The ONLY holder of the full MAC key. In real SPDZ the full alpha never
// exists: each party samples its own alpha_p locally (alpha = Σ alpha_p is
// never reconstructed), and a MAC share of alpha*x is assembled from the
// pairwise cross products x_i*alpha_j, obtained via correlated OT (MASCOT) or
// homomorphic encryption (Overdrive), plus the local alpha_p*x_p term.
// Replacing this struct with such an offline is the SPDZ de-cheat: protocol
// code only ever touches dealer.my_alpha and dealer.deal(v), both of which a
// real offline can serve without the `alpha` field existing.
template <int nP> struct FakeDealer {
  int party;
  std::mt19937 rng;  // shared seed: every party runs the identical stream
  uint32_t alpha;    // CHEAT: full key, derived from the shared seed
  uint32_t my_alpha; // this party's key share — all that protocols may see

  FakeDealer(int party, uint32_t seed) : party(party), rng(seed) {
    alpha = rng() % Q;
    my_alpha = arith_split<nP>(rng, alpha)[party];
  }
  // This party's authenticated share of a dealer-chosen value v.
  AuthShare deal(uint32_t v) { return auth_split<nP>(rng, v, alpha)[party]; }
};

} // namespace mldsa
#endif // MLDSA_SPDZ_H
