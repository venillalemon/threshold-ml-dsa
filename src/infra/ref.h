// src/infra/ref.h — ML-DSA parameters, F_q / R_q utilities, and the plaintext
// reference algorithms (FIPS 204), shared by every layer.
//
// Everything here is PUBLIC LOCAL computation — no shares, no circuits, no
// network. Three of the algorithms are CHEAT stubs of FIPS 204's SHAKE-based
// ones (expand_a, h256_stub, sample_in_ball_stub); swapping in real Keccak
// changes no protocol message besides the values themselves.
#ifndef MLDSA_REF_H
#define MLDSA_REF_H

#include <array>
#include <cstdint>
#include <random>
#include <vector>

namespace mldsa {

// ---- parameter sets ---------------------------------------------------------
enum ParamSet { MLDSA_44, MLDSA_65, MLDSA_87 };

#ifndef MLDSA_PARAM
#define MLDSA_PARAM MLDSA_65
#endif
constexpr ParamSet PARAM = MLDSA_PARAM;

constexpr double lambda = 32.0; // statistical parameter for oversampled batches
constexpr int32_t Q = 8380417;
constexpr int32_t L = 23; // bits of a coefficient in [0, Q)
constexpr int N = 256;    // R_q = Z_q[x]/(x^N+1)
constexpr int K = PARAM == MLDSA_44 ? 4 : (PARAM == MLDSA_65 ? 6 : 8);
constexpr int ELL = PARAM == MLDSA_44 ? 4 : (PARAM == MLDSA_65 ? 5 : 7);
constexpr int COEFF_COUNT = N * K;       // a length-k vector of polynomials
constexpr int Y_COEFF_COUNT = N * ELL;   // a length-ell vector of polynomials
constexpr int32_t G2 = (PARAM == MLDSA_44) ? (Q - 1) / 88 : (Q - 1) / 32; // gamma2
constexpr int ETA = PARAM == MLDSA_65 ? 4 : 2;
constexpr int TAU = PARAM == MLDSA_44 ? 39 : (PARAM == MLDSA_65 ? 49 : 60);   // +-1s in c
constexpr int OMEGA = PARAM == MLDSA_44 ? 80 : (PARAM == MLDSA_65 ? 55 : 75); // max hint weight
constexpr int BETA = TAU * ETA; // ||c*s||_inf, ||c*e||_inf <= BETA
constexpr int GAMMA1 = 1 << 19;
constexpr int Y_WIDTH = 20;      // log2(GAMMA1) + 1 for 2-complement
constexpr int POW2ROUND_D = 13;  // Power2Round dropped bits d

inline int param_k(ParamSet p) {
  return p == MLDSA_44 ? 4 : (p == MLDSA_65 ? 6 : 8);
}
inline int param_coeffs(ParamSet p) {
  return 256 * param_k(p);
}
inline const char* param_name(ParamSet p) {
  return p == MLDSA_44 ? "ML-DSA-44" : (p == MLDSA_65 ? "ML-DSA-65" : "ML-DSA-87");
}

// ---- F_q scalar ops ---------------------------------------------------------
inline uint32_t fq_add(uint32_t a, uint32_t b) {
  return (uint32_t)(((uint64_t)a + b) % Q);
}
inline uint32_t fq_sub(uint32_t a, uint32_t b) {
  return (uint32_t)(((uint64_t)a + Q - b) % Q);
}
inline uint32_t fq_mul(uint32_t a, uint32_t b) {
  return (uint32_t)(((uint64_t)a * b) % Q);
}
inline uint32_t modq_i64(int64_t v) {
  const int64_t m = v % Q;
  return (uint32_t)(m < 0 ? m + Q : m);
}

// ---- public R_q arithmetic --------------------------------------------------

// acc[N] += c * x over the integers in Z[X]/(X^N+1); c sparse in {-1,0,1}.
inline void cpoly_mul_acc(const std::vector<int32_t>& c, const int64_t* x, int64_t* acc) {
  for (int j = 0; j < N; ++j) {
    if (!c[j])
      continue;
    for (int d = 0; d < N; ++d) {
      int deg = j + d;
      int64_t s = (int64_t)c[j] * x[d];
      if (deg >= N) {
        deg -= N;
        s = -s;
      }
      acc[deg] += s;
    }
  }
}

// public A (rows x cols x N polys) times public x (cols*N in [0,Q)), mod q.
inline std::vector<uint32_t> matvec_pub_negacyclic(const std::vector<uint32_t>& A,
                                                   const std::vector<uint32_t>& x, int rows,
                                                   int cols) {
  std::vector<uint32_t> out((size_t)rows * N, 0);
  for (int row = 0; row < rows; ++row)
    for (int col = 0; col < cols; ++col)
      for (int a_deg = 0; a_deg < N; ++a_deg) {
        const uint32_t a = A[((size_t)row * cols + col) * N + a_deg];
        if (!a)
          continue;
        for (int x_deg = 0; x_deg < N; ++x_deg) {
          const int degree = a_deg + x_deg;
          const int out_deg = degree < N ? degree : degree - N;
          const uint32_t a_mod_q = degree < N ? a : (uint32_t)Q - a;
          out[row * N + out_deg] =
              fq_add(out[row * N + out_deg], fq_mul(a_mod_q, x[col * N + x_deg]));
        }
      }
  return out;
}

// ---- FIPS 204 reference algorithms (plaintext) ------------------------------

// Decompose: r = r1*(2*g2) + r0 with r0 in (-g2, g2], plus the q-1 special case.
inline void ref_decompose(int32_t r, int32_t g2, int32_t& r1, int32_t& r0) {
  const int32_t a = 2 * g2;
  int32_t rp = r % a;
  r0 = (rp > g2) ? rp - a : rp;
  if (r - r0 == Q - 1) {
    r1 = 0;
    r0 -= 1;
  } else
    r1 = (r - r0) / a;
}
inline int32_t ref_highbits(int32_t r, int32_t g2) {
  int32_t r1, r0;
  ref_decompose(r, g2, r1, r0);
  return r1;
}
inline int32_t ref_lowbits(int32_t r, int32_t g2) {
  int32_t r1, r0;
  ref_decompose(r, g2, r1, r0);
  return r0;
}

// Power2Round: t = t1*2^d + t0 over the integers for t in [0,Q),
// with t0 in (-2^(d-1), 2^(d-1)].
inline void ref_power2round(uint32_t t, int32_t& t1, int32_t& t0) {
  t1 = ((int32_t)t + (1 << (POW2ROUND_D - 1)) - 1) >> POW2ROUND_D;
  t0 = (int32_t)t - (t1 << POW2ROUND_D);
}

// MakeHint (signer side): does adding z change the high bits of r?
// h = [HighBits(r) != HighBits(r+z)]; z may be negative (e.g. -c*t0).
inline int ref_makehint(int64_t z, int32_t r, int32_t g2) {
  return ref_highbits(r, g2) != ref_highbits((int32_t)modq_i64((int64_t)r + z), g2);
}

// UseHint (verifier side).
inline int32_t ref_usehint(int h, int32_t r, int32_t g2) {
  const int32_t m = (Q - 1) / (2 * g2);
  int32_t r1, r0;
  ref_decompose(r, g2, r1, r0);
  if (!h)
    return r1;
  return r0 > 0 ? (r1 + 1) % m : (r1 - 1 + m) % m;
}

// ---- SHAKE-shaped stubs (CHEATS) -------------------------------------------

// CHEAT(Fcoin): all parties must agree on fresh public randomness; a fixed
// constant keeps every party consistent without a coin-flip protocol.
inline std::array<uint32_t, 8> fcoin_rho() {
  return {0x9e3779b9u, 0x7f4a7c15u, 0x85ebca6bu, 0xc2b2ae35u,
          0x27d4eb2fu, 0x165667b1u, 0xd3a2646cu, 0xfd7046c5u};
}

// CHEAT(ExpandA): rho-seeded mt19937 instead of SHAKE-128. Deterministic in
// rho, so every party derives the same A locally — the property KeyGen/Sign
// actually rely on. Rejection to [0,Q) matches FIPS coefficient sampling.
inline void expand_a(const std::array<uint32_t, 8>& rho, std::vector<uint32_t>& A) {
  std::seed_seq seq(rho.begin(), rho.end());
  std::mt19937 gen(seq);
  A.resize((size_t)K * ELL * N);
  for (uint32_t& a : A)
    a = gen() % Q;
}

// CHEAT(H): 256-bit hash via seed_seq mixing; real ML-DSA uses SHAKE-256.
inline void h256_stub(const std::vector<uint32_t>& input, std::array<uint32_t, 8>& out) {
  std::seed_seq seq(input.begin(), input.end());
  seq.generate(out.begin(), out.end());
}

// CHEAT(H): SampleInBall driven by a seeded mt19937 instead of a SHAKE-256
// stream. The Fisher-Yates structure (TAU signed swaps) is FIPS 204's.
inline void sample_in_ball_stub(const std::array<uint32_t, 8>& seed, std::vector<int32_t>& c) {
  std::seed_seq seq(seed.begin(), seed.end());
  std::mt19937 gen(seq);
  c.assign(N, 0);
  for (int i = N - TAU; i < N; ++i) {
    const int j = gen() % (i + 1);
    const int s = gen() & 1;
    c[i] = c[j];
    c[j] = s ? -1 : 1;
  }
}

} // namespace mldsa
#endif // MLDSA_REF_H
