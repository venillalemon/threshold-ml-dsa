// src/protocol/prepsign.h — Pi_PrepSign, offline (pre-challenge).
//
// Arithmetic (edaBits, w = A*y + e_w, the field opening of delta_w = w + R2) is
// BDOZ over F_q via the FakeDealer, exactly as before. The Boolean half — the
// A2B recovery of w and its Decompose — is now ONE compiled circuit evaluated
// under authenticated GMW (the offline circuit domain, C_ad below), producing
// public w1 and retained authenticated shares of w0 that feed the boundary
// producers in sign.h / sign_2round.h.
#ifndef MLDSA_PREPSIGN_H
#define MLDSA_PREPSIGN_H

#include "a2b.h"        // sub_modq kernel
#include "decompose.h"  // decompose<PARAM>, OW0, OW1
#include "dealer.h"
#include "phase1_util.h"

#include <cstdint>
#include <cstdio>
#include <vector>

namespace mldsa {

#ifdef TEST
inline std::vector<uint32_t> test_opened_y, test_opened_w;
#endif

// C_ad: per coefficient, w = sub_modq(c, R2) then (w1, w0) = Decompose(w).
// Inputs  [R2 : L*COEFF_COUNT][c : L*COEFF_COUNT]  (c enters as public bits).
// Outputs [w1 : OW1*COEFF_COUNT][w0 : OW0*COEFF_COUNT].
inline const emp::circuit::BooleanProgram& recover_decompose_program() {
  static const emp::circuit::BooleanProgram prog = [] {
    emp::RecordCtx ctx;
    const uint32_t base = ctx.external_input((uint32_t)(2 * L * COEFF_COUNT));
    const uint32_t c_base = base + (uint32_t)(L * COEFF_COUNT);
    std::vector<emp::UInt_T<emp::RecordCtx, OW1>> w1s;
    std::vector<emp::UInt_T<emp::RecordCtx, OW0>> w0s;
    w1s.reserve(COEFF_COUNT);
    w0s.reserve(COEFF_COUNT);
    for (int i = 0; i < COEFF_COUNT; ++i) {
      auto r2 = in_uint<L>(ctx, base + (uint32_t)(L * i));
      auto c = in_uint<L>(ctx, c_base + (uint32_t)(L * i));
      auto w = sub_modq(ctx, c, r2);
      auto d = decompose<PARAM>(ctx, w);
      w1s.push_back(d.w1);
      w0s.push_back(d.w0);
    }
    std::vector<emp::RecordCtx::Wire> outw;
    outw.reserve((size_t)(OW1 + OW0) * COEFF_COUNT);
    for (auto& v : w1s) {
      emp::RecordCtx::Wire w[OW1];
      v.pack_wires(w);
      outw.insert(outw.end(), w, w + OW1);
    }
    for (auto& v : w0s) {
      emp::RecordCtx::Wire w[OW0];
      v.pack_wires(w);
      outw.insert(outw.end(), w, w + OW0);
    }
    ctx.finish(outw);
    return std::move(ctx.prog);
  }();
  return prog;
}

template <int nP> struct PrepSignOut {
  emp::wrk::ShareVec<nP> w0_2;         // OW0-bit GMW shares per coefficient
  SharePair<nP, Y_WIDTH, true> y;      // Boolean = bin(R_y), arithmetic = y = gamma1 - R_y
  std::vector<uint32_t> w1;            // public HighBits
};

template <int nP>
inline PrepSignOut<nP> prepsign(Backend<nP>& bk, int party, FakeDealer<nP>& dealer,
                                const std::vector<uint32_t>& A) {
  PrepSignOut<nP> out;

  // <y>_2,<y>_q and <e_w>_q. <e_w>_2 is intentionally unused.
  out.y = dealer.template deal_y_edabit<Y_WIDTH>(Y_COEFF_COUNT);
  auto ew = dealer.deal_secret_poly(COEFF_COUNT, ETA); // e_w: arithmetic only

#ifdef TEST
  auto open_q = [&](const std::vector<FqShare<nP>>& s) {
    return open_fq_checked(bk.io(), party, dealer.my_alpha_f, s);
  };
  auto centered = [](uint32_t x) { return x > (uint32_t)Q / 2 ? (int32_t)x - Q : (int32_t)x; };
  const std::vector<uint32_t> opened_y = open_q(out.y.fq_share);
  const std::vector<uint32_t> opened_ew = open_q(ew.fq);
  int y_bad = 0, ew_bad = 0;
  for (uint32_t x : opened_y)
    y_bad += centered(x) <= -GAMMA1 || centered(x) > GAMMA1;
  for (uint32_t x : opened_ew)
    ew_bad += centered(x) < -ETA || centered(x) > ETA;
  emp::expecting(y_bad == 0, "PrepSign test: y out of (-GAMMA1,GAMMA1]");
  emp::expecting(ew_bad == 0, "PrepSign test: e_w out of [-ETA,ETA]");
#endif

  // <w>_q = A<y>_q + <e_w>_q.
  std::vector<FqShare<nP>> w_share = ew.fq;
  matvec_negacyclic(A, out.y.fq_share, w_share, K, ELL);

  // A2B mask R2, and the checked field open c = w + R2.
  SharePair<nP, L, true> r2 = dealer.template deal_fq_edabit<L>(COEFF_COUNT);
  std::vector<FqShare<nP>> c_share((size_t)COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i)
    c_share[(size_t)i] = w_share[(size_t)i] + r2.fq_share[(size_t)i];
#ifdef TAMPER_C
  if (party == 1 && COEFF_COUNT > 0)
    c_share[0].val ^= 1;
#endif
  const std::vector<uint32_t> c = open_fq_checked(bk.io(), party, dealer.my_alpha_f, c_share);

  // C_ad inputs: [R2 boolean bits][c public bits].
  emp::wrk::ShareVec<nP> in;
  in.reserve((size_t)2 * L * COEFF_COUNT);
  in.insert(in.end(), r2.two_share.begin(), r2.two_share.end());
  for (int i = 0; i < COEFF_COUNT; ++i)
    push_public_word(bk, in, c[(size_t)i], L);

  auto outs = bk.gmw().evaluate(recover_decompose_program(), in);
  emp::expecting((int)outs.size() == (OW1 + OW0) * COEFF_COUNT, "prepsign: C_ad output width");

  // Open w1; retain w0 shares.
  emp::wrk::ShareVec<nP> w1shares(outs.begin(), outs.begin() + (size_t)OW1 * COEFF_COUNT);
  const std::vector<uint8_t> w1bits = bk.gmw().open(w1shares, "prepsign-w1");
  out.w1.assign((size_t)COEFF_COUNT, 0);
  for (int i = 0; i < COEFF_COUNT; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < OW1; ++k)
      v |= (uint32_t)w1bits[(size_t)(i * OW1 + k)] << k;
    out.w1[(size_t)i] = v;
  }
  out.w0_2.assign(outs.begin() + (size_t)OW1 * COEFF_COUNT, outs.end());

#ifdef TEST
  const std::vector<uint32_t> opened_w = open_q(w_share);
  emp::wrk::ShareVec<nP> w0shares(outs.begin() + (size_t)OW1 * COEFF_COUNT, outs.end());
  const std::vector<uint8_t> w0bits = bk.gmw().open(w0shares, "prepsign-w0-test");
  int bad = 0;
  for (int i = 0; i < COEFF_COUNT; ++i) {
    int32_t w1e, w0e;
    ref_decompose((int32_t)opened_w[(size_t)i], G2, w1e, w0e);
    int32_t w0v = 0;
    for (int k = 0; k < OW0; ++k)
      w0v |= (int32_t)w0bits[(size_t)(i * OW0 + k)] << k;
    w0v -= (w0v >> (OW0 - 1)) << OW0; // sign-extend two's complement
    if ((int32_t)out.w1[(size_t)i] != w1e || w0v != w0e)
      ++bad;
  }
  if (party == 1)
    std::printf("%s  nP=%d  %d coefficients  ->  %d wrong\n", param_name(PARAM), nP, COEFF_COUNT, bad);
  emp::expecting(bad == 0, "PrepSign test: w0/w1 differ from FIPS Decompose");
  test_opened_y = opened_y;
  test_opened_w = opened_w;
#endif

  return out;
}

} // namespace mldsa

#endif // MLDSA_PREPSIGN_H
