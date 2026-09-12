// src/protocol/sign_2round.h — Pi_MLDSA Sign, single-shot (T = 1), late-challenge
// window design WITHOUT slots (mldsa44_late_challenge_functionalities.pdf,
// Figs. 4-5). This is the baseline the slot-restricted variant in sign.h is
// compared against; build it with -DMLDSA_SLOT=OFF (make SLOT=0).
//
// The challenge-dependent displacement is tiny: with H^z = beta - c*s and
// H^w = beta + c*e, 0 <= H <= 2*beta, and for U^z = gamma1 - y, U^w = gamma2 - w0
//   accept_z  <=>  2*beta < U^z + H^z < 2*gamma1        (= strict ||z|| bound)
//   accept_w  <=>  2*beta < U^w + H^w < 2*gamma2        (= strict ||r0|| bound)
// so the decision only depends on whether U is within 2*beta of a boundary
// (known BEFORE c) and on a RING_K-bit window of U + H. Everything wide is
// therefore pre-challenge:
//
//   Phase 1 (offline, ordinary Boolean chunk):  per coefficient, from Boolean
//     U and the Boolean half r of a ring eDaBit (R^, r): lower = [U <= 2beta],
//     upper = [U >= span - 2beta], edge = lower ^ upper,
//     B = base(U) (Eq. 4), J = shift(U, span) (Eq. 3), both selected on the
//     widened power-of-two edge regions -- see window_producer below
//     M = R + J mod 2^RING_K
//   Offline:   sess.prepare(Phase 2), a fixed circuit on (B, M^z, M^w, p).
//   Online:    H^ = beta -/+ c*<s/e>_ring (LOCAL, linear on ring shares),
//              open p^ = H^ - R^ over Z_{2^RING_W}, BDOZ-checked        (1 flight)
//              run(prepared): T = M + p mod 2^RING_K, accept_i = T[RING_K-1],
//              Q = B + low(T) mod 2^Y_WIDTH (= gamma1 + beta - z),
//              r = AND_i accept_i, T^z = r ? Q : gamma1+beta;   reveal (r, T^z)
//   Local:     z = gamma1 + beta - T^z, then MakeHint etc.
//
// Correctness of the window (all mod 2^RING_K, HALF = 2^(RING_K-1), 2beta < HALF):
//   lower:  T = U + H + HALF - 2beta - 1, top bit set  <=>  U + H >= 2beta + 1
//   upper:  T = U + H - span, top bit set              <=>  U + H < span
//   middle: T = HALF + H, always accepted (U itself is safely inside)
// and in each accepted case B + low(T) = U + H exactly.
#ifndef MLDSA_SIGN_2ROUND_H
#define MLDSA_SIGN_2ROUND_H

#include "circuit.h"
#include "edabits.h"
#include "keygen.h"
#include "prepsign.h"
#include "ref.h"
#include "spdz.h"
#include <emp-ag/emp-ag.h>
#include <emp-tool/circuits/frontend/frontend.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace mldsa {

constexpr int HALF = 1 << (RING_K - 1);
constexpr int64_t SPAN_Z = 2 * (int64_t)GAMMA1;    // U^z in [0, SPAN_Z)
constexpr int64_t SPAN_W = 2 * (int64_t)G2;        // U^w in [0, SPAN_W)
constexpr int64_t C_LOW = 2 * BETA + 1;            // B in the lower case
constexpr int64_t C_J = HALF - 2 * BETA - 1;       // J offset in the lower case
constexpr int64_t B_UP = SPAN_Z - HALF;            // B in the upper case
constexpr int64_t Z_OFF = (int64_t)GAMMA1 + BETA;  // T^z = Z_OFF - z
static_assert(SPAN_Z == (int64_t{1} << Y_WIDTH), "U^z must be exactly Y_WIDTH bits");
static_assert(SPAN_W < (int64_t{1} << OW0), "U^w must fit OW0 bits");
static_assert(SPAN_Z % (1 << RING_K) == 0, "SPAN_Z must be a multiple of the window modulus");

constexpr int BW = Y_COEFF_COUNT * Y_WIDTH;                // B, all response coefficients
constexpr int MZW = Y_COEFF_COUNT * RING_K;                // M^z
constexpr int MWW = COEFF_COUNT * RING_K;                  // M^w
constexpr int PW = (Y_COEFF_COUNT + COEFF_COUNT) * RING_K; // public p, z first then w
constexpr int OUTW = 1 + BW;                               // r || T^z

template <class Ctx> using BV = emp::BitVec_T<Ctx, BW>;
template <class Ctx> using MZV = emp::BitVec_T<Ctx, MZW>;
template <class Ctx> using MWV = emp::BitVec_T<Ctx, MWW>;
template <class Ctx> using PV = emp::BitVec_T<Ctx, PW>;
template <class Ctx> using OUTV = emp::BitVec_T<Ctx, OUTW>;

// (c, z, h) with explicit bot markers: z_ok=false => (c,bot,bot) [r0 or z
// rejected]; h_ok=false => (c,z,bot). All public — this IS the signature.
struct Signature {
  std::vector<int32_t> c; // challenge polynomial, N coefficients in {-1,0,1}
  bool z_ok = false;
  std::vector<int32_t> z; // N*ELL coefficients, |z| <= GAMMA1-BETA (iff z_ok)
  bool h_ok = false;
  std::vector<uint8_t> h; // N*K hint bits (iff h_ok)
};

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

// ---- Phase 1: boundary producer (offline Boolean) ---------------------------

// Specialised producer (paper Sec. 4.1). L = 2^b, D = 2L, ny = b+1 = RING_K.
//
// The edge regions are widened to the power-of-two boundary, which is legal
// because shift()'s lower formula is valid for every U <= L (then
// L-2beta-1 <= J+H <= L+2beta-1 < D and K[b]=1 iff U+H >= 2beta+1) and its
// upper formula for every U >= S-L; the interior L <= U < S-L is accepted for
// all reachable H. So
//     lower = [U < L]        = the top UW-b bits are all zero      (NOR tree)
//     upper = [U >= S - L]   = the top UW-b bits are all one       (AND tree,
//                              when S = 2^UW; otherwise a compare)
// Widening also makes base()'s high UW-b bits equal U's in all three cases, so
// only the low b bits of B are selected, and the byte X = edge & low_b(U) is
// shared between B and J:
//     low_b(B) = low_b(U) ^ X ^ (2beta+1) * lower      (free given X)
//     low_b(J) = X + (L - 2beta - 1) * lower  (mod 2^b)
//     J[b]     = !lower                     (S = 2^UW; else edge ? U[b] : 1)
//     M        = R + J                      (mod D)
// Cost: 2(UW-b-1) + b + (b-1) + b = 2*UW + b - 3 AND for a power-of-two span
// (41 for MLDSA44's response row, matching the paper), one more for the J[b]
// of a non-power-of-two span, and a UW-bit compare instead of the AND tree.
template <int UW, int64_t SPAN, class Ctx>
inline void window_producer(Ctx& ctx, const Bit_T<Ctx>* U, const Bit_T<Ctx>* R,
                            Bit_T<Ctx>* B, Bit_T<Ctx>* M) {
  using Bit = Bit_T<Ctx>;
  constexpr int b = RING_K - 1;                    // L = 2^b
  constexpr bool SPAN_POW2 = SPAN == (int64_t{1} << UW);
  constexpr int64_t S0 = SPAN & (int64_t)RING_KMASK; // SPAN mod D
  static_assert((S0 & ((int64_t{1} << b) - 1)) == 0,
                "SPAN mod D must have zero low bits so low_b(J) needs no correction");
  static_assert(2 * BETA < (int64_t{1} << b), "L must exceed 2*beta");

  Bit hi[UW - b];
  for (int k = 0; k < UW - b; ++k)
    hi[k] = U[b + k];
  const Bit lower = !or_tree(hi, UW - b);
  Bit upper = Bit::constant(ctx, false);
  if constexpr (SPAN_POW2)
    upper = and_tree(hi, UW - b);
  else
    upper = ge_const_u<UW>(ctx, U, (uint64_t)(SPAN - (int64_t{1} << b)));
  const Bit edge = lower ^ upper;

  Bit X[b];
  for (int k = 0; k < b; ++k)
    X[k] = edge & U[k];

  if (B != nullptr) {
    for (int k = 0; k < b; ++k)
      B[k] = ((C_LOW >> k) & 1) ? (U[k] ^ X[k] ^ lower) : (U[k] ^ X[k]);
    for (int k = b; k < UW; ++k)
      B[k] = U[k]; // base()'s high bits coincide with U's in every region
  }

  // low_b(J) = X + C_J * lower, with its carry folded into bit b. The widened
  // lower region reaches U = L - 1, so U + C_J can exceed 2^b - 1 and the carry
  // is real (with the paper's exact [U <= 2beta] flag it never is).
  Bit J[RING_K], carry = Bit::constant(ctx, false);
  for (int k = 0; k < b; ++k) {
    if (((C_J >> k) & 1) != 0)
      fa(J[k], carry, X[k], lower, carry);
    else
      ha(J[k], carry, X[k], carry);
  }
  if constexpr (SPAN_POW2)
    J[b] = !lower ^ carry; // interior and upper both leave U[b] = 1
  else
    J[b] = !(edge & !U[b]) ^ carry; // (edge ? U[b] : 1) + carry
  if constexpr ((S0 >> b) & 1)
    J[b] = J[b] ^ upper; // subtracting S0 = L flips bit b (MLDSA65 low-bit row)
  add_ripple<RING_K>(ctx, M, R, J, false);
}

// ---- Phase 2: the fixed post-challenge circuit -------------------------------

template <class Ctx>
inline OUTV<Ctx> phase2_body(Ctx& ctx, const BV<Ctx>& B, const MZV<Ctx>& Mz, const MWV<Ctx>& Mw,
                             const PV<Ctx>& P) {
  using Bit = Bit_T<Ctx>;
  std::vector<Bit> ok((size_t)Y_COEFF_COUNT + COEFF_COUNT);
  std::vector<Bit> Q((size_t)BW);
  Bit m[RING_K], p[RING_K], T[RING_K];

  for (int idx = 0; idx < Y_COEFF_COUNT; ++idx) {
    for (int k = 0; k < RING_K; ++k) {
      m[k] = Mz[idx * RING_K + k];
      p[k] = P[idx * RING_K + k];
    }
    add_ripple<RING_K>(ctx, T, m, p, false);
    ok[(size_t)idx] = T[RING_K - 1];
    Bit b[Y_WIDTH], ext[Y_WIDTH];
    for (int k = 0; k < Y_WIDTH; ++k) {
      b[k] = B[idx * Y_WIDTH + k];
      ext[k] = k < RING_K - 1 ? T[k] : Bit::constant(ctx, false); // low(T), zero-extended
    }
    add_ripple<Y_WIDTH>(ctx, &Q[(size_t)idx * Y_WIDTH], b, ext, false);
  }
  for (int idx = 0; idx < COEFF_COUNT; ++idx) {
    for (int k = 0; k < RING_K; ++k) {
      m[k] = Mw[idx * RING_K + k];
      p[k] = P[(Y_COEFF_COUNT + idx) * RING_K + k];
    }
    add_ripple<RING_K>(ctx, T, m, p, false);
    ok[(size_t)(Y_COEFF_COUNT + idx)] = T[RING_K - 1];
  }

  const Bit r = and_tree(ok.data(), (int)ok.size());
  std::vector<Bit> outb((size_t)OUTW);
  outb[0] = r;
  for (int t = 0; t < BW; ++t)
    outb[(size_t)(1 + t)] = r & Q[(size_t)t]; // Z_out = r * ((B + low_b(K^z)) mod 2^a)
  return OUTV<Ctx>::from_bit_values(ctx, outb.data());
}

// Compiled once per process: the circuit depends on the parameter set only.
inline const auto& phase2_circuit() {
  using emp::RecordCtx;
  static const auto c =
      emp::frontend::compile<BV<RecordCtx>, MZV<RecordCtx>, MWV<RecordCtx>, PV<RecordCtx>>(
          [](RecordCtx& ctx, BV<RecordCtx> B, MZV<RecordCtx> Mz, MWV<RecordCtx> Mw,
             PV<RecordCtx> P) { return phase2_body(ctx, B, Mz, Mw, P); });
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

  // Ring eDaBits (R^, r) for every response and low-bit position.
  RingEdabits<nP> rz = ring_edabits<nP, RING_K>(sess, party, dealer, Y_COEFF_COUNT);
  RingEdabits<nP> rw = ring_edabits<nP, RING_K>(sess, party, dealer, COEFF_COUNT);

  // Phase 1: Boolean boundary producers -> B, M^z, M^w (stay authenticated).
  std::vector<Wire> bw((size_t)BW), mzw((size_t)MZW), mww((size_t)MWW);
  {
    Bit U[Y_WIDTH], R[RING_K], B[Y_WIDTH], M[RING_K];
    for (int idx = 0; idx < Y_COEFF_COUNT; ++idx) {
      // <y>_2 stores v = y - 1 (two's complement); U = gamma1 - y = (gamma1-1) - v
      // = ~v with the top bit restored -- pure rewiring.
      const auto v = sess.template adopt_authenticated_input<emp::UInt_T<Ctx, Y_WIDTH>>(
          &y.two_share[(size_t)Y_WIDTH * idx]);
      for (int k = 0; k < Y_WIDTH; ++k)
        U[k] = k < Y_WIDTH - 1 ? !v[k] : v[k];
      const auto r = sess.template adopt_authenticated_input<UR>(&rz.two_share[(size_t)RING_K * idx]);
      for (int k = 0; k < RING_K; ++k)
        R[k] = r[k];
      window_producer<Y_WIDTH, SPAN_Z>(sess.ctx(), U, R, B, M);
      for (int k = 0; k < Y_WIDTH; ++k)
        B[k].pack_wires(&bw[(size_t)idx * Y_WIDTH + k]);
      for (int k = 0; k < RING_K; ++k)
        M[k].pack_wires(&mzw[(size_t)idx * RING_K + k]);
    }
  }
  {
    Bit w0b[OW0], nw0[OW0], U[OW0], R[RING_K], M[RING_K];
    for (int idx = 0; idx < COEFF_COUNT; ++idx) {
      unpack_bits<OW0>(w0_2[(size_t)idx], w0b);
      for (int k = 0; k < OW0; ++k)
        nw0[k] = !w0b[k];
      add_const_ripple<OW0>(sess.ctx(), U, nw0, (uint64_t)G2, true); // U = gamma2 - w0
      const auto r = sess.template adopt_authenticated_input<UR>(&rw.two_share[(size_t)RING_K * idx]);
      for (int k = 0; k < RING_K; ++k)
        R[k] = r[k];
      window_producer<OW0, SPAN_W>(sess.ctx(), U, R, static_cast<Bit*>(nullptr), M);
      for (int k = 0; k < RING_K; ++k)
        M[k].pack_wires(&mww[(size_t)idx * RING_K + k]);
    }
  }
  const BV<Ctx> B_v = BV<Ctx>::from_wires(sess.ctx(), bw.data());
  const MZV<Ctx> Mz_v = MZV<Ctx>::from_wires(sess.ctx(), mzw.data());
  const MWV<Ctx> Mw_v = MWV<Ctx>::from_wires(sess.ctx(), mww.data());
  sess.checkpoint(); // run Phase 1 now; B / M become soldering sources
  if (party == 1) // G1 split: prepsign (edaBits + A2B + Decompose) vs producers
    std::printf("G1SPLIT prepsign=%llu producers=%llu\n", (unsigned long long)(and_a1 - and_a0),
                (unsigned long long)(sess.num_and() - and_a1));
#ifdef TEST
  timer_lap("ring edabits + phase 1");
#endif

  // Garble the fixed Phase 2 circuit before the message exists.
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

  // H^z = beta - c*s, H^w = beta + c*e on the ring shares (local, linear);
  // p^ = H^ - R^ over Z_{2^RING_W}, one checked BDOZ open for all positions.
  std::vector<RingShare<nP>> ph((size_t)Y_COEFF_COUNT + COEFF_COUNT);
  {
    std::vector<RingShare<nP>> cs((size_t)Y_COEFF_COUNT), ce((size_t)COEFF_COUNT);
    for (int p = 0; p < ELL; ++p)
      cpoly_mul_ring(c_nz, &kp.s.ring_share[(size_t)p * N], &cs[(size_t)p * N]);
    for (int p = 0; p < K; ++p)
      cpoly_mul_ring(c_nz, &kp.e.ring_share[(size_t)p * N], &ce[(size_t)p * N]);
    for (int i = 0; i < Y_COEFF_COUNT; ++i)
      ph[(size_t)i] = add_public_ring(RingShare<nP>{} - cs[(size_t)i], (uint64_t)BETA, party,
                                      dealer.my_alpha_r) -
                      rz.ring_share[(size_t)i];
    for (int i = 0; i < COEFF_COUNT; ++i)
      ph[(size_t)Y_COEFF_COUNT + i] =
          add_public_ring(ce[(size_t)i], (uint64_t)BETA, party, dealer.my_alpha_r) -
          rw.ring_share[(size_t)i];
  }
#ifdef TAMPER_C
  if (party == 1)
    ph[0].val.w[0] ^= 1;
#endif
  const std::vector<RingVal> phat = open_ring_checked(sess.io(), party, dealer.my_alpha_r, ph);

  // p = p^ mod 2^RING_K as public circuit inputs.
  typename PV<Ctx>::clear_t pb{};
  for (size_t i = 0; i < phat.size(); ++i)
    for (int k = 0; k < RING_K; ++k)
      pb[i * RING_K + (size_t)k] = ((phat[i].w[0] >> k) & 1) != 0;
  const PV<Ctx> P_v = sess.template input<PV<Ctx>>(emp::PUBLIC, pb);

  // One prepared run, one batched reveal.
  const OUTV<Ctx> out = sess.run(std::move(prepared), B_v, Mz_v, Mw_v, P_v);
  const auto opened = sess.reveal(out, 1); // Open((r, Z_out), P1): only the evaluator decodes
  if (online_rounds_out)
    *online_rounds_out = 3; // flight 1, run, reveal
  if (party != 1)
    return sig;
  emp::expecting(opened.has_value(), "Sign: online reveal failed");
#ifdef TEST
  timer_lap("open + run + reveal");
#endif

  const auto& ob = opened.value();
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
      emp::expecting(!ob[(size_t)t], "Sign test: rejected run released a nonzero Z_out bit");
#endif

  if (!accept)
    return sig; // (c, bot, bot)

  sig.z.resize(Y_COEFF_COUNT);
  for (int i = 0; i < Y_COEFF_COUNT; ++i) {
    int64_t t = 0; // Z_out,i = gamma1 + beta - z_i
    for (int k = 0; k < Y_WIDTH; ++k)
      t |= (int64_t)ob[(size_t)(1 + i * Y_WIDTH + k)] << k;
    sig.z[i] = (int32_t)(Z_OFF - t);
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
#endif // MLDSA_SIGN_2ROUND_H
