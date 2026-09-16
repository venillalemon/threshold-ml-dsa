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

#include "backend.h"
#include "circuit.h"
#include "dealer.h"
#include "keygen.h"
#include "phase1_util.h"
#include "a2b.h"
#include "decompose.h"
#include "ref.h"
#include "wrk_phase2.h"
#include <emp-tool/circuits/frontend/frontend.h>
#include <emp-tool/ir/context/record.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace mldsa {

constexpr int M_TOTAL = Y_COEFF_COUNT + COEFF_COUNT; // response + low-bit positions
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

// The compiled C_post program (fixed inputs B||Mz||Mw, late public P = delta_H).
inline const emp::circuit::BooleanProgram& phase2_program() { return phase2_circuit().program(); }

// C_pre port layout (both modes): secret [R_y][R2][R_H], then public [c].
constexpr uint32_t CP_RY = 0;
constexpr uint32_t CP_R2 = CP_RY + (uint32_t)(Y_WIDTH * Y_COEFF_COUNT);
constexpr uint32_t CP_RH = CP_R2 + (uint32_t)(L * COEFF_COUNT);
constexpr uint32_t CP_C = CP_RH + (uint32_t)(RING_K * M_TOTAL);
constexpr uint32_t CP_IN = CP_C + (uint32_t)(L * COEFF_COUNT);
constexpr int W1W = OW1 * COEFF_COUNT; // opened HighBits
constexpr size_t CP_W0_OFF = (size_t)W1W + BW + MZW + MWW; // TEST-only w0 output offset

// C_pre (paper step 3), ONE offline GMW circuit: recover w = c - R2, Decompose,
// then the boundary producers.
// Inputs  [R_y : Y_WIDTH*nz][R2 : L*nw][R_H^z : RING_K*nz][R_H^w : RING_K*nw][c : L*nw public]
// Outputs [w1 : W1W (opened)][B : BW][Mz : MZW][Mw : MWW]   (+ [w0 : OW0*nw] under TEST)
inline const emp::circuit::BooleanProgram& cpre_program() {
  using Ctx = emp::RecordCtx;
  using Bit = emp::Bit_T<Ctx>;
  using Wire = Ctx::Wire;
  static const emp::circuit::BooleanProgram prog = [] {
    Ctx ctx;
    const uint32_t base = ctx.external_input(CP_IN);
    std::vector<Wire> outW1, outB, outMz, outMw, outW0;
    outW1.reserve(W1W);
    outB.reserve(BW);
    outMz.reserve(MZW);
    outMw.reserve(MWW);
    Bit U[Y_WIDTH > OW0 ? Y_WIDTH : OW0], R[RING_K], B[Y_WIDTH], M[RING_K];
    for (int i = 0; i < Y_COEFF_COUNT; ++i) { // response rows: U = R_y directly
      for (int k = 0; k < Y_WIDTH; ++k)
        U[k] = wire_bit(ctx, base + CP_RY + (uint32_t)(Y_WIDTH * i + k));
      for (int k = 0; k < RING_K; ++k)
        R[k] = wire_bit(ctx, base + CP_RH + (uint32_t)(RING_K * i + k));
      window_producer<Y_WIDTH, SPAN_Z>(ctx, U, R, B, M);
      for (int k = 0; k < Y_WIDTH; ++k) {
        Wire w;
        B[k].pack_wires(&w);
        outB.push_back(w);
      }
      for (int k = 0; k < RING_K; ++k) {
        Wire w;
        M[k].pack_wires(&w);
        outMz.push_back(w);
      }
    }
    for (int i = 0; i < COEFF_COUNT; ++i) { // low-bit rows
        auto r2 = in_uint<L>(ctx, base + CP_R2 + (uint32_t)(L * i));
        auto cc = in_uint<L>(ctx, base + CP_C + (uint32_t)(L * i));
        auto d = decompose<PARAM>(ctx, sub_modq(ctx, cc, r2)); // paper step 3: w, then (w1, w0)
        {
          Wire ww[OW1];
          d.w1.pack_wires(ww);
          outW1.insert(outW1.end(), ww, ww + OW1);
        }
        Bit w0b[OW0], nw0[OW0];
        unpack_bits<OW0>(d.w0, w0b);
#ifdef TEST
        {
          Wire ww[OW0];
          d.w0.pack_wires(ww);
          outW0.insert(outW0.end(), ww, ww + OW0);
        }
#endif
        for (int k = 0; k < OW0; ++k)
          nw0[k] = !w0b[k];
        add_const_ripple<OW0>(ctx, U, nw0, (uint64_t)G2, true); // U = gamma2 - w0
        for (int k = 0; k < RING_K; ++k)
          R[k] = wire_bit(ctx, base + CP_RH + (uint32_t)(RING_K * (Y_COEFF_COUNT + i) + k));
        window_producer<OW0, SPAN_W>(ctx, U, R, static_cast<Bit*>(nullptr), M);
        for (int k = 0; k < RING_K; ++k) {
          Wire w;
          M[k].pack_wires(&w);
          outMw.push_back(w);
        }
    }
    std::vector<Wire> outs;
    outs.insert(outs.end(), outW1.begin(), outW1.end());
    outs.insert(outs.end(), outB.begin(), outB.end());
    outs.insert(outs.end(), outMz.begin(), outMz.end());
    outs.insert(outs.end(), outMw.begin(), outMw.end());
    outs.insert(outs.end(), outW0.begin(), outW0.end());
    ctx.finish(outs);
    return std::move(ctx.prog);
  }();
  return prog;
}

// No-op default for sign()'s after_prepsign hook (see below).
struct NoopHook {
  void operator()() const {}
};

// ---- Sign -------------------------------------------------------------------
// `after_prepsign` fires at the offline/online boundary (after the eDaBits, the
// eDaBits, Phase 1 and the prepared garbling). `online_rounds_out` receives the
// number of online synchronization barriers; `g2_ands_out` the AND count of the
// prepared post-challenge circuit.
template <int nP, class Hook = NoopHook>
inline Signature sign(Backend<nP>& bk, int party, FakeDealer<nP>& dealer, KeyPair<nP>& kp,
                      const std::vector<uint32_t>& msg, Hook after_prepsign = Hook{},
                      int64_t* online_rounds_out = nullptr, int64_t* g2_ands_out = nullptr,
                      int64_t* g1_ands_out = nullptr) {
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
  // Paper step 1: the eDaBit families (dealer-dealt).
  SharePair<nP, Y_WIDTH, true> y = dealer.template deal_y_edabit<Y_WIDTH>(Y_COEFF_COUNT);
  auto ew = dealer.deal_secret_poly(COEFF_COUNT, ETA); // e_w: arithmetic only
  SharePair<nP, L, true> r2 = dealer.template deal_fq_edabit<L>(COEFF_COUNT); // R_w
  RingEdabits<nP> rz = dealer.template deal_ring_edabit<RING_K>(Y_COEFF_COUNT); // R_H^z
  RingEdabits<nP> rw = dealer.template deal_ring_edabit<RING_K>(COEFF_COUNT);   // R_H^w

  // Paper step 2: w = A y + e_w on F_q shares; open delta_w as c = w + R2
  // (one checked BDOZ flight; the circuit recovers w = c - R2).
  std::vector<FqShare<nP>> w_share = ew.fq;
  matvec_negacyclic(kp.A, y.fq_share, w_share, K, ELL);
  std::vector<FqShare<nP>> c_share((size_t)COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i)
    c_share[(size_t)i] = w_share[(size_t)i] + r2.fq_share[(size_t)i];
#ifdef TAMPER_C
  if (party == 1)
    c_share[0].val ^= 1;
#endif
  const std::vector<uint32_t> c = open_fq_checked(bk.io(), party, dealer.my_alpha_f, c_share);
#ifdef TEST
  auto open_q = [&](const std::vector<FqShare<nP>>& s) {
    return open_fq_checked(bk.io(), party, dealer.my_alpha_f, s);
  };
  auto centered = [](uint32_t x) { return x > (uint32_t)Q / 2 ? (int32_t)x - Q : (int32_t)x; };
  std::vector<uint32_t> test_opened_y = open_q(y.fq_share), test_opened_w = open_q(w_share);
  {
    const std::vector<uint32_t> opened_ew = open_q(ew.fq);
    int y_bad = 0, ew_bad = 0;
    for (uint32_t x : test_opened_y)
      y_bad += centered(x) <= -GAMMA1 || centered(x) > GAMMA1;
    for (uint32_t x : opened_ew)
      ew_bad += centered(x) < -ETA || centered(x) > ETA;
    emp::expecting(y_bad == 0, "Sign test: y out of (-GAMMA1,GAMMA1]");
    emp::expecting(ew_bad == 0, "Sign test: e_w out of [-ETA,ETA]");
  }
  timer_lap("edabits + w + open c");
#endif

  // Paper step 3: ONE offline GMW circuit C_pre.
  const uint64_t and_a0 = bk.gmw().ands_evaluated;
  emp::wrk::ShareVec<nP> pin;
  pin.reserve((size_t)CP_IN);
  pin.insert(pin.end(), y.two_share.begin(), y.two_share.end());
  pin.insert(pin.end(), r2.two_share.begin(), r2.two_share.end());
  pin.insert(pin.end(), rz.two_share.begin(), rz.two_share.end());
  pin.insert(pin.end(), rw.two_share.begin(), rw.two_share.end());
  for (int i = 0; i < COEFF_COUNT; ++i)
    push_public_word(bk, pin, c[(size_t)i], L);
  emp::expecting(pin.size() == (size_t)CP_IN, "Sign: C_pre input width");
  emp::wrk::ShareVec<nP> outs = bk.gmw().evaluate(cpre_program(), pin);
  pin.clear();
  pin.shrink_to_fit();
  emp::expecting(outs.size() >= CP_W0_OFF, "Sign: C_pre output width");
  // Open w1 (needed for the challenge); B, Mz, Mw stay authenticated shares.
  std::vector<uint32_t> w1((size_t)COEFF_COUNT, 0);
  {
    emp::wrk::ShareVec<nP> w1s(outs.begin(), outs.begin() + W1W);
    const std::vector<uint8_t> w1bits = bk.gmw().open(w1s, "cpre-w1");
    for (int i = 0; i < COEFF_COUNT; ++i)
      for (int k = 0; k < OW1; ++k)
        w1[(size_t)i] |= (uint32_t)w1bits[(size_t)(i * OW1 + k)] << k;
  }
  emp::wrk::ShareVec<nP> bm(outs.begin() + W1W, outs.begin() + W1W + (size_t)(BW + MZW + MWW));
  if (g1_ands_out)
    *g1_ands_out = (int64_t)(bk.gmw().ands_evaluated - and_a0);

#ifdef TEST
  { // Decompose oracle: the circuit's (w1, w0) must equal FIPS Decompose(w).
    emp::wrk::ShareVec<nP> w0s(outs.begin() + CP_W0_OFF, outs.begin() + CP_W0_OFF + (size_t)OW0 * COEFF_COUNT);
    const std::vector<uint8_t> w0bits = bk.gmw().open(w0s, "cpre-w0-test");
    int bad = 0;
    for (int i = 0; i < COEFF_COUNT; ++i) {
      int32_t w1e, w0e;
      ref_decompose((int32_t)test_opened_w[(size_t)i], G2, w1e, w0e);
      int32_t w0v = 0;
      for (int k = 0; k < OW0; ++k)
        w0v |= (int32_t)w0bits[(size_t)(i * OW0 + k)] << k;
      w0v -= (w0v >> (OW0 - 1)) << OW0;
      if ((int32_t)w1[(size_t)i] != w1e || w0v != w0e)
        ++bad;
    }
    if (party == 1)
      std::printf("%s  nP=%d  %d coefficients  ->  %d wrong\n", param_name(PARAM), nP, COEFF_COUNT, bad);
    emp::expecting(bad == 0, "Sign test: w0/w1 differ from FIPS Decompose");
  }
  timer_lap("C_pre (GMW)");
#endif

  // Garble C_post OFFLINE: rows, fixed-input labels and output masks all
  // delivered to P1 before the message exists (wrk_phase2.h).
  if (g2_ands_out)
    *g2_ands_out = (int64_t)emp::wrk::count_gc_ands(phase2_program());
  WrkOffline<nP> off = wrk_offline<nP>(bk, phase2_program(), bm, (uint32_t)(BW + MZW + MWW));
  outs.clear();
  outs.shrink_to_fit();
#ifdef TEST
  timer_lap("wrk offline garble");
#endif
  bk.gmw().finish(); // complete all malicious-COT checks; no GMW past this point
  after_prepsign();  // offline/online boundary -- everything below needs msg

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

  // Flight 1: H^z = beta - c*s, H^w = beta + c*e on the ring shares (local);
  // p^ = H^ - R^ over Z_{2^RING_W}, one checked BDOZ open for all positions.
  std::vector<RingShare<nP>> ph((size_t)M_TOTAL);
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
  const std::vector<RingVal> phat = open_ring_checked(bk.io(), party, dealer.my_alpha_r, ph);

  // Late public bits P = p^ mod 2^RING_K, z block then w block.
  std::vector<uint8_t> late((size_t)PW);
  for (size_t i = 0; i < phat.size(); ++i)
    for (int k = 0; k < RING_K; ++k)
      late[i * RING_K + (size_t)k] = (uint8_t)((phat[i].w[0] >> k) & 1);

  // Flight 2: garblers send selected labels; P1 evaluates and decodes locally.
  auto res = wrk_online<nP>(bk, off, late);
  if (online_rounds_out)
    *online_rounds_out = 2; // flight 1 (delta_H open), flight 2 (labels)
  if (party != 1)
    return sig;
  emp::expecting(res.has_value(), "Sign: online evaluation failed");
#ifdef TEST
  timer_lap("flight1 open + flight2 eval");
#endif

  const std::vector<uint8_t>& ob = *res;
  const bool accept = ob[0];

#ifdef TEST
  // Plaintext oracle: recompute r0 / z from the opened e, s (keygen) and w, y
  // (C_pre); the circuit's accept bit and (if accepted) z must agree.
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
