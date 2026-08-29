// src/protocol/prepsign.h — Pi_PrepSign orchestration.
#ifndef MLDSA_PREPSIGN_H
#define MLDSA_PREPSIGN_H

#include "a2b.h"
#include "decompose.h"
#include "edabits.h"
#include "rand.h"
#include "spdz.h"
#include <emp-ag/emp-ag.h>

#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace mldsa {

#ifdef TEST
// TEST-only oracle: plaintext values opened by the latest prepsign() call, so
// downstream TEST checks (sign.h) can cross-check circuits against plaintext.
inline std::vector<uint32_t> test_opened_y, test_opened_w;
#endif

// PrepSign outputs. w0_2 remains a vector of secret AG circuit values so it can
// be consumed directly by the next circuit. y_2 is the authenticated Boolean
// sharing returned by y_edabits (the 20-bit two's-complement encoding of y-1).
// w1 is public at every party.
template <int nP> using PrepSignCtx = typename emp::AGMPCSession<nP>::ctx_t;

// A is the public matrix ExpandA(rho) from KeyGen (keygen.h), flattened
// K x ELL x N. The paper has Sign recompute it locally from rho; here the
// caller passes KeyPair::A.
template <int nP>
inline void prepsign(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer,
                     const std::vector<uint32_t>& A,
                     std::vector<emp::UInt_T<PrepSignCtx<nP>, OW0>>& w0_2,
                     emp::ag::AShareBundleVec<nP>& y_2, std::vector<uint32_t>& w1) {
  using Ctx = PrepSignCtx<nP>;
  using U23 = emp::UInt_T<Ctx, L>;

#ifdef TEST
  auto lap_start = emp::clock_start();
  auto timer_lap = [&](const char* label) {
    if (party == 1)
      std::printf("  [timer] %-20s %8.1f ms\n", label, emp::time_from(lap_start) / 1000.0);
    lap_start = emp::clock_start();
  };
#endif

  // <y>_2,<y>_q and <e_w>_q. <e_w>_2 is intentionally unused.
  SharePair<nP, Y_WIDTH, false> y = y_edabits<nP, Y_WIDTH>(sess, party, dealer, Y_COEFF_COUNT);
#ifdef TEST
  timer_lap("y edabits");
#endif
  SharePair<nP, (ETA == 2 ? 3 : 4), false> ew =
      rand_edabits<nP, ETA, COEFF_COUNT>(sess, party, dealer);
#ifdef TEST
  timer_lap("e_w edabits");
#endif

#ifdef TEST
  // Test-only range checks; every opened arithmetic vector is MAC-checked.
  auto open_q_shares = [&](const std::vector<AuthShare>& shares) {
    std::vector<uint32_t> val(shares.size()), mac(shares.size());
    for (size_t i = 0; i < shares.size(); ++i) {
      val[i] = shares[i].val;
      mac[i] = shares[i].mac;
    }
    std::vector<uint32_t> opened = open_additive_modq(sess.io(), party, val);
    spdz_maccheck(sess.io(), party, dealer.my_alpha, mac, opened);
    return opened;
  };
  auto centered = [](uint32_t x) { return x > (uint32_t)Q / 2 ? (int32_t)x - Q : (int32_t)x; };

  const std::vector<uint32_t> opened_y = open_q_shares(y.q_share);
  int y_min = GAMMA1 + 1, y_max = -GAMMA1 - 1, y_bad = 0;
  for (uint32_t x : opened_y) {
    const int value = centered(x);
    y_min = value < y_min ? value : y_min;
    y_max = value > y_max ? value : y_max;
    y_bad += value <= -GAMMA1 || value > GAMMA1;
  }

  const std::vector<uint32_t> opened_ew = open_q_shares(ew.q_share);
  int ew_min = ETA + 1, ew_max = -ETA - 1, ew_bad = 0;
  for (uint32_t x : opened_ew) {
    const int value = centered(x);
    ew_min = value < ew_min ? value : ew_min;
    ew_max = value > ew_max ? value : ew_max;
    ew_bad += value < -ETA || value > ETA;
  }
  if (party == 1) {
    std::printf("  y range [%d, %d], %d outside (-%d, %d]\n", y_min, y_max, y_bad, GAMMA1, GAMMA1);
    std::printf("  e_w range [%d, %d], %d outside [-%d, %d]\n", ew_min, ew_max, ew_bad, ETA, ETA);
  }
  emp::expecting(y_bad == 0, "PrepSign test: opened y coefficient outside (-GAMMA1,GAMMA1]");
  emp::expecting(ew_bad == 0, "PrepSign test: opened e_w coefficient outside [-ETA,ETA]");
  timer_lap("test range checks");
#endif

  // <w>_q=A<y>_q+<e_w>_q in F_q[x]/(x^N+1).
  std::vector<AuthShare> w_share = ew.q_share;
  matvec_negacyclic(A, y.q_share, w_share, K, ELL); // add A*y onto e_w
#ifdef TEST
  timer_lap("A*y + e_w (local)");
#endif

  // <w>_q -> <w>_2 via F_A2B (a2b.h: mask with an edaBit r, open c = w + r
  // under a MACCheck, recompute w = c - r inside the circuit).
  std::vector<U23> w_2 = a2b<nP>(sess, party, dealer, w_share);
#ifdef TEST
  timer_lap("A2B");
#endif

  // Decompose using AG
  std::vector<emp::UInt_T<Ctx, OW1>> w1_wires;
  w1_wires.reserve(COEFF_COUNT);
  std::vector<emp::UInt_T<Ctx, OW0>> w0_wires;
  w0_wires.reserve(COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i) {
    auto d = decompose<PARAM>(sess.ctx(), w_2[i]);
    w1_wires.push_back(d.w1);
    w0_wires.push_back(d.w0);
  }
#ifdef TEST
  timer_lap("decompose circuit");
#endif

  w1.resize(COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i) {
    const auto public_w1 = sess.reveal(w1_wires[i], emp::PUBLIC);
    emp::expecting(public_w1.has_value(), "PrepSign: public w1 reveal failed");
    w1[i] = (uint32_t)public_w1.value();
  }
#ifdef TEST
  timer_lap("w1 reveal");
#endif

#ifdef TEST
  // Test oracle input. This opening does not exist in the protocol build.
  std::vector<uint32_t> w_val(COEFF_COUNT), w_mac(COEFF_COUNT);
  std::vector<int32_t> opened_w0(COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i) {
    w_val[i] = w_share[i].val;
    w_mac[i] = w_share[i].mac;
    const auto public_w0 = sess.reveal(w0_wires[i], emp::PUBLIC);
    emp::expecting(public_w0.has_value(), "PrepSign test: public w0 reveal failed");
    const int32_t value = (int32_t)public_w0.value();
    opened_w0[i] = value - ((value >> (OW0 - 1)) << OW0);
  }
  const std::vector<uint32_t> opened_w = open_additive_modq(sess.io(), party, w_val);
  spdz_maccheck(sess.io(), party, dealer.my_alpha, w_mac, opened_w);
  int bad = 0;
  for (int i = 0; i < COEFF_COUNT; ++i) {
    int32_t expected_w1, expected_w0;
    ref_decompose((int32_t)opened_w[i], G2, expected_w1, expected_w0);
    if ((int32_t)w1[i] != expected_w1 || opened_w0[i] != expected_w0)
      ++bad;
  }
  if (party == 1)
    std::printf("%s  nP=%d  %d coefficients  ->  %d wrong\n", param_name(PARAM), nP,
                COEFF_COUNT, bad);
  emp::expecting(bad == 0, "PrepSign test: w0/w1 differ from FIPS Decompose");
  test_opened_y = opened_y;
  test_opened_w = opened_w;
  timer_lap("test w0/w1 oracle");
#endif

  w0_2 = std::move(w0_wires);
  y_2 = std::move(y.two_share);
}

} // namespace mldsa

#endif // MLDSA_PREPSIGN_H
