// src/protocol/sign_slot.h — Pi_MLDSA Sign, single-shot (T = 1), SLOT-RESTRICTED
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
#ifndef MLDSA_SIGN_SLOT_H
#define MLDSA_SIGN_SLOT_H

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
constexpr int DHW = NS * RING_K;                  // routed residues d^
constexpr int RHOW = Y_COEFF_COUNT * L;           // bin23(rho), the F_q response pad
constexpr int UW_ALL = NS * M_TOTAL;              // routing matrix u[j][i]
constexpr int OUTW = 1 + RHOW;                    // r || Z_out
// C_post port layout: fixed [M^][bin23(rho)][u], then routed late [d^].
constexpr uint32_t P2_MH = 0;
constexpr uint32_t P2_RHO = P2_MH + MHW;
constexpr uint32_t P2_U = P2_RHO + RHOW;
constexpr uint32_t P2_DH = P2_U + UW_ALL;
constexpr uint32_t P2_FIXED = P2_DH;              // == MHW + RHOW + UW_ALL
constexpr uint32_t P2_IN = P2_DH + DHW;

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

// ---- Phase 2: the post-challenge circuit C_post -------------------------------
//
// Inputs  [M^ : MHW][bin23(rho) : RHOW][u : NS*M_TOTAL]  (fixed, installed offline)
//         [d^ : NS*RING_K]                              (routed late, see wrk_online_routed)
// Outputs [r][Z_out : RHOW]
//
// u is never read by a gate: it is installed so P1 holds its labels and every
// party its masks, and d^[j][l] = XOR_{i: delta_H,i[l]=1} u[j][i] is formed by
// free XOR online (paper: "the evaluator forms it by XOR-ing active labels it
// already holds"). The d^ ports carry fresh offline masks; flight 2 is the
// re-masking that moves P1 from the XORed u wires onto those ports (one
// 1+16n-byte block per port per garbler, wrk_online_routed). Every gate
// (K_j = M^_j + d^_j, r = AND_j K_j[top], Z_out = r * bin23(rho)) is garbled
// offline.
inline const emp::circuit::BooleanProgram& phase2_program() {
  using Ctx = emp::RecordCtx;
  using Bit = emp::Bit_T<Ctx>;
  using Wire = Ctx::Wire;
  static const emp::circuit::BooleanProgram prog = [] {
    Ctx ctx;
    const uint32_t base = ctx.external_input((uint32_t)P2_IN);
    // Slot adders + conjunction.
    std::vector<Bit> ok((size_t)NS);
    Bit a[RING_K], b[RING_K], Kk[RING_K];
    for (int j = 0; j < NS; ++j) {
      for (int k = 0; k < RING_K; ++k) {
        a[k] = wire_bit(ctx, base + P2_MH + (uint32_t)(j * RING_K + k));
        b[k] = wire_bit(ctx, base + P2_DH + (uint32_t)(j * RING_K + k));
      }
      add_ripple<RING_K>(ctx, Kk, a, b, false);
      ok[(size_t)j] = Kk[RING_K - 1];
    }
    const Bit r = and_tree(ok.data(), NS);
    std::vector<Wire> outs;
    outs.reserve((size_t)OUTW);
    {
      Wire w;
      r.pack_wires(&w);
      outs.push_back(w);
    }
    for (int t = 0; t < RHOW; ++t) {
      const Bit z = r & wire_bit(ctx, base + P2_RHO + (uint32_t)t); // Z_out = r * bin(rho)
      Wire w;
      z.pack_wires(&w);
      outs.push_back(w);
    }
    ctx.finish(outs);
    return std::move(ctx.prog);
  }();
  return prog;
}

// C_pre port layout (both modes): secret [R_y][R2][R_H], then public [c].
constexpr uint32_t CP_RY = 0;
constexpr uint32_t CP_R2 = CP_RY + (uint32_t)(Y_WIDTH * Y_COEFF_COUNT);
constexpr uint32_t CP_RH = CP_R2 + (uint32_t)(L * COEFF_COUNT);
constexpr uint32_t CP_C = CP_RH + (uint32_t)(RING_K * M_TOTAL);
constexpr uint32_t CP_IN = CP_C + (uint32_t)(L * COEFF_COUNT);
constexpr int W1W = OW1 * COEFF_COUNT; // opened HighBits
constexpr size_t CP_W0_OFF = (size_t)W1W + MHW + UW_ALL + 1; // TEST-only w0 output offset

// C_pre (slot, paper step 3), ONE offline GMW circuit: recover w = c - R2,
// Decompose, producers -> (M_i, edge_i), running one-hot decoder -> u, shift
// prefix, order-preserving compaction of M into the slots, empty-slot MUX, ovf.
// Inputs  [R_y : Y_WIDTH*nz][R2 : L*nw][R_H : RING_K*M_TOTAL][c : L*nw public]
// Outputs [w1 : W1W (opened)][M^ : MHW][u : NS*M_TOTAL, u[j][i] at j*M_TOTAL+i][ovf : 1 (opened)]
//         (+ [w0 : OW0*nw] under TEST)
// u is retained: online, d^_j[l] = XOR_{i : delta_H,i[l]=1} u[j][i] is a free XOR
// of its labels with PUBLIC coefficients (paper Fig. 1 routing).
inline const emp::circuit::BooleanProgram& cpre_program() {
  using Ctx = emp::RecordCtx;
  using Bit = emp::Bit_T<Ctx>;
  using Wire = Ctx::Wire;
  static const emp::circuit::BooleanProgram prog = [] {
    Ctx ctx;
    const uint32_t base = ctx.external_input(CP_IN);
    std::vector<Wire> outW1, outW0;
    outW1.reserve(W1W);

    // Producers -> (M_i, edge_i).
    std::vector<std::array<Bit, RING_K>> Marr((size_t)M_TOTAL);
    std::vector<std::array<Bit, DW>> Sarr((size_t)M_TOTAL);
    std::vector<Bit> edge((size_t)M_TOTAL);
    {
      Bit U[Y_WIDTH > OW0 ? Y_WIDTH : OW0], R[RING_K];
      for (int idx = 0; idx < M_TOTAL; ++idx) {
        for (int k = 0; k < RING_K; ++k)
          R[k] = wire_bit(ctx, base + CP_RH + (uint32_t)(RING_K * idx + k));
        if (idx < Y_COEFF_COUNT) {
          for (int k = 0; k < Y_WIDTH; ++k) // the eDaBit's Boolean half IS R_y = U^z
            U[k] = wire_bit(ctx, base + CP_RY + (uint32_t)(Y_WIDTH * idx + k));
          window_producer<Y_WIDTH>(ctx, U, R, SPAN_Z, edge[(size_t)idx], Marr[(size_t)idx].data());
        } else {
          const int i = idx - Y_COEFF_COUNT;
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
          window_producer<OW0>(ctx, U, R, SPAN_W, edge[(size_t)idx], Marr[(size_t)idx].data());
        }
      }
    }

    // One-hot decoder: u[j][i] = edge_i AND [cnt_i = j]; shift prefix; ovf.
    std::vector<std::vector<Bit>> u((size_t)NS, std::vector<Bit>((size_t)M_TOTAL));
    Bit ovf;
    {
      std::vector<Bit> h((size_t)NS + 1, Bit::constant(ctx, false));
      h[0] = Bit::constant(ctx, true);
      std::vector<Bit> ucol((size_t)NS + 1);
      Bit shift[DW];
      for (int k = 0; k < DW; ++k)
        shift[k] = Bit::constant(ctx, false);
      for (int i = 0; i < M_TOTAL; ++i) {
        for (int j = 0; j <= NS; ++j)
          ucol[(size_t)j] = edge[(size_t)i] & h[(size_t)j];
        for (int j = 0; j < NS; ++j)
          u[(size_t)j][(size_t)i] = ucol[(size_t)j];
        for (int j = NS; j >= 0; --j)
          h[(size_t)j] = h[(size_t)j] ^ (j > 0 ? ucol[(size_t)j - 1] : Bit::constant(ctx, false)) ^
                         ucol[(size_t)j];
        for (int k = 0; k < DW; ++k)
          Sarr[(size_t)i][(size_t)k] = edge[(size_t)i] & shift[k];
        add_bit_ripple<DW>(ctx, shift, shift, !edge[(size_t)i]);
      }
      Bit acc = h[0];
      for (int j = 1; j <= NS; ++j)
        acc = acc ^ h[(size_t)j];
      ovf = !acc;
    }

    compact<RING_K>(Marr, Sarr, M_TOTAL); // M moves into the slots

    // Outputs: w1, then M^ with the empty-slot MUX to L = 2^(RING_K-1), then u, ovf.
    std::vector<Wire> outs;
    outs.reserve((size_t)W1W + (size_t)MHW + (size_t)NS * M_TOTAL + 1 + outW0.size());
    outs.insert(outs.end(), outW1.begin(), outW1.end());
    for (int j = 0; j < NS; ++j) {
      Bit empty = Bit::constant(ctx, true);
      for (int i = 0; i < M_TOTAL; ++i)
        empty = empty ^ u[(size_t)j][(size_t)i];
      for (int k = 0; k < RING_K; ++k) {
        const Bit m = Marr[(size_t)j][(size_t)k];
        const Bit bit = k == RING_K - 1 ? (m ^ (empty & !m)) : (m & !empty);
        Wire w;
        bit.pack_wires(&w);
        outs.push_back(w);
      }
    }
    for (int j = 0; j < NS; ++j)
      for (int i = 0; i < M_TOTAL; ++i) {
        Wire w;
        u[(size_t)j][(size_t)i].pack_wires(&w);
        outs.push_back(w);
      }
    {
      Wire w;
      ovf.pack_wires(&w);
      outs.push_back(w);
    }
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
  RingEdabits<nP> rh = dealer.template deal_ring_edabit<RING_K>(M_TOTAL);          // R_H
  SharePair<nP, L, true> rho = dealer.template deal_fq_edabit<L>(Y_COEFF_COUNT);    // response pad

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
  pin.insert(pin.end(), rh.two_share.begin(), rh.two_share.end());
  for (int i = 0; i < COEFF_COUNT; ++i)
    push_public_word(bk, pin, c[(size_t)i], L);
  emp::expecting(pin.size() == (size_t)CP_IN, "Sign: C_pre input width");
  emp::wrk::ShareVec<nP> outs = bk.gmw().evaluate(cpre_program(), pin);
  pin.clear();
  pin.shrink_to_fit();
  emp::expecting(outs.size() >= CP_W0_OFF, "Sign: C_pre output width");
  // Paper: Open((w1, ovf)) together. Retire the attempt on overflow.
  std::vector<uint32_t> w1((size_t)COEFF_COUNT, 0);
  {
    emp::wrk::ShareVec<nP> both(outs.begin(), outs.begin() + W1W);
    both.push_back(outs[(size_t)W1W + MHW + UW_ALL]); // ovf
    const std::vector<uint8_t> bits = bk.gmw().open(both, "cpre-w1-ovf");
    for (int i = 0; i < COEFF_COUNT; ++i)
      for (int k = 0; k < OW1; ++k)
        w1[(size_t)i] |= (uint32_t)bits[(size_t)(i * OW1 + k)] << k;
    if (bits[(size_t)W1W]) {
      std::fprintf(stderr, "Sign: slot overflow (>%d edge coefficients), retire attempt\n", NS);
      bk.gmw().finish();
      return sig; // (c, bot, bot) -- happens with probability ~2^-207
    }
  }
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

  // Garble C_post OFFLINE. Fixed inputs = M^ || bin23(rho) || u (all installed
  // now); d^ is a ROUTED late input formed online by free XOR of u.
  emp::wrk::ShareVec<nP> fixed(outs.begin() + W1W, outs.begin() + W1W + MHW);
  fixed.insert(fixed.end(), rho.two_share.begin(), rho.two_share.end());
  fixed.insert(fixed.end(), outs.begin() + W1W + MHW, outs.begin() + W1W + MHW + (size_t)UW_ALL);
  outs.clear();
  outs.shrink_to_fit();
  if (g2_ands_out)
    *g2_ands_out = (int64_t)emp::wrk::count_gc_ands(phase2_program());
  WrkOffline<nP> off = wrk_offline<nP>(bk, phase2_program(), fixed, P2_FIXED, /*late_routed=*/true);
  fixed.clear();
  fixed.shrink_to_fit();
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

  // Flight 1 (paper step 5): one opening carrying both domains --
  //   [[delta_H]]_D = (beta*1 - c[[s]]_D, beta*1 + c[[e]]_D) - [[R_H]]_D   (ring)
  //   [[V]]_q       = (gamma1+beta)*1 - [[y]]_q - c[[s]]_q + [[rho]]_q     (field)
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
  open_ring_and_fq(bk.io(), party, dealer.my_alpha_r, dealer.my_alpha_f, dh, vbar_sh, opened, vbar);

  // Routing: S[l] = { i : delta_H,i[l] = 1 } -- the PUBLIC coefficients of
  // d^_j[l] = XOR_{i in S[l]} u[j][i]. Flight 2: garblers send the re-masking
  // blocks for the d^ ports; P1 forms d^ by free XOR, evaluates, decodes.
  std::vector<std::vector<int>> S((size_t)RING_K);
  for (int l = 0; l < RING_K; ++l)
    for (int i = 0; i < M_TOTAL; ++i)
      if ((opened[(size_t)i].w[0] >> l) & 1)
        S[(size_t)l].push_back(i);
  auto res = wrk_online_routed<nP>(bk, off, P2_U, (uint32_t)M_TOTAL, P2_DH, (uint32_t)RING_K,
                                   (uint32_t)NS, S);
  if (online_rounds_out)
    *online_rounds_out = 2; // flight 1 (delta_H, V), flight 2 (d^ re-masking)
  if (party != 1)
    return sig;
  emp::expecting(res.has_value(), "Sign: online evaluation failed");
#ifdef TEST
  timer_lap("flight1 + routing + eval");
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
#endif // MLDSA_SIGN_SLOT_H
