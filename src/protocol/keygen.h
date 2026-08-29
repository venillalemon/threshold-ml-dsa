// src/protocol/keygen.h — Pi_MLDSA KeyGen.
//
// Paper steps -> code:
//   1. rho <- Fcoin; <s>,<e> <- F_smallnormpoly   fcoin_rho (CHEAT), rand_edabits
//   2. A = ExpandA(rho); <t>_q = A<s>_q + <e>_q   expand_a (CHEAT), matvec_negacyclic
//      Fopen(<t>_q)                               open_additive_modq + spdz_maccheck
//   3. (t1,t0) = Power2Round(t,d); tr = H(rho,t1) ref_power2round, h256_stub (CHEAT)
//   4. store (rho,tr,t,t1,t0,<s>_2,<e>_2), pk=(rho,t1)   KeyPair
//
// The Fcoin / ExpandA / H stand-ins (fcoin_rho, expand_a, h256_stub) and
// ref_power2round live in ref.h — all public local computations.
#ifndef MLDSA_KEYGEN_H
#define MLDSA_KEYGEN_H

#include "rand.h"
#include "ref.h"
#include "spdz.h"
#include <emp-ag/emp-ag.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace mldsa {

// Everything KeyGen stores. pk = (rho, t1); tr seeds mu in Sign. The paper
// stores only the Boolean halves <s>_2, <e>_2 — the q-shares inside the
// SharePairs are dead after t is opened, kept only because SharePair owns both.
template <int nP> struct KeyPair {
  std::array<uint32_t, 8> rho;
  std::array<uint32_t, 8> tr;
  std::vector<uint32_t> A;         // ExpandA(rho), cached (Sign recomputes it in the paper)
  std::vector<uint32_t> t;         // opened t, coefficients in [0,Q)
  std::vector<int32_t> t1, t0;     // Power2Round(t, d)
  SharePair<nP, (ETA == 2 ? 3 : 4), false> s; // <s>_2 (N*ELL coefficients)
  SharePair<nP, (ETA == 2 ? 3 : 4), false> e; // <e>_2 (N*K coefficients)
#ifdef TEST
  std::vector<uint32_t> opened_s, opened_e; // TEST-only plaintext oracle
#endif
};

template <int nP>
inline KeyPair<nP> keygen(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer) {
  KeyPair<nP> kp;

#ifdef TEST
  auto lap_start = emp::clock_start();
  auto timer_lap = [&](const char* label) {
    if (party == 1)
      std::printf("  [timer] %-20s %8.1f ms\n", label, emp::time_from(lap_start) / 1000.0);
    lap_start = emp::clock_start();
  };
#endif

  // 1. rho from Fcoin; A = ExpandA(rho) locally at every party.
  kp.rho = fcoin_rho();
  expand_a(kp.rho, kp.A);

  // <s>: ell polynomials, <e>: k polynomials, both eta-small (F_smallnormpoly).
  kp.s = rand_edabits<nP, ETA, Y_COEFF_COUNT>(sess, party, dealer);
#ifdef TEST
  timer_lap("s smallnormpoly");
#endif
  kp.e = rand_edabits<nP, ETA, COEFF_COUNT>(sess, party, dealer);
#ifdef TEST
  timer_lap("e smallnormpoly");
#endif

  // 2. <t>_q = A<s>_q + <e>_q, then Fopen(<t>_q) gated by a MACCheck.
  std::vector<AuthShare> t_share = kp.e.q_share;
  matvec_negacyclic(kp.A, kp.s.q_share, t_share, K, ELL);
#ifdef TEST
  timer_lap("A*s + e (local)");
#endif

  std::vector<uint32_t> t_val(COEFF_COUNT), t_mac(COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i) {
    t_val[i] = t_share[i].val;
    t_mac[i] = t_share[i].mac;
  }
  kp.t = open_additive_modq(sess.io(), party, t_val);
  spdz_maccheck(sess.io(), party, dealer.my_alpha, t_mac, kp.t);
#ifdef TEST
  timer_lap("open t + MACCheck");
#endif

  // 3. Power2Round and tr; both local on the now-public t.
  kp.t1.resize(COEFF_COUNT);
  kp.t0.resize(COEFF_COUNT);
  for (int i = 0; i < COEFF_COUNT; ++i)
    ref_power2round(kp.t[i], kp.t1[i], kp.t0[i]);

  std::vector<uint32_t> tr_in(kp.rho.begin(), kp.rho.end());
  for (int32_t v : kp.t1)
    tr_in.push_back((uint32_t)v);
  h256_stub(tr_in, kp.tr);

#ifdef TEST
  // Test-only: open s,e (destroys secrecy — TEST builds only) and check
  // range, t = A*s + e in the clear, and the Power2Round identity.
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

  const std::vector<uint32_t> opened_s = open_q_shares(kp.s.q_share);
  const std::vector<uint32_t> opened_e = open_q_shares(kp.e.q_share);
  kp.opened_s = opened_s;
  kp.opened_e = opened_e;
  int s_bad = 0, e_bad = 0;
  for (uint32_t x : opened_s)
    s_bad += centered(x) < -ETA || centered(x) > ETA;
  for (uint32_t x : opened_e)
    e_bad += centered(x) < -ETA || centered(x) > ETA;
  emp::expecting(s_bad == 0, "KeyGen test: opened s coefficient outside [-ETA,ETA]");
  emp::expecting(e_bad == 0, "KeyGen test: opened e coefficient outside [-ETA,ETA]");

  // Plaintext A*s + e must reproduce the opened t.
  std::vector<uint32_t> t_ref = matvec_pub_negacyclic(kp.A, opened_s, K, ELL);
  for (int i = 0; i < COEFF_COUNT; ++i)
    t_ref[i] = fq_add(t_ref[i], opened_e[i]);
  int t_bad = 0, p2r_bad = 0;
  for (int i = 0; i < COEFF_COUNT; ++i) {
    t_bad += t_ref[i] != kp.t[i];
    const bool identity = (int64_t)kp.t1[i] * (1 << POW2ROUND_D) + kp.t0[i] == (int32_t)kp.t[i];
    const bool range = kp.t0[i] > -(1 << (POW2ROUND_D - 1)) && kp.t0[i] <= (1 << (POW2ROUND_D - 1));
    p2r_bad += !identity || !range;
  }
  if (party == 1)
    std::printf("  KeyGen: t mismatches %d, Power2Round violations %d\n", t_bad, p2r_bad);
  emp::expecting(t_bad == 0, "KeyGen test: opened t != A*s + e");
  emp::expecting(p2r_bad == 0, "KeyGen test: Power2Round identity/range violated");
  timer_lap("test keygen checks");
#endif

  return kp;
}

} // namespace mldsa
#endif // MLDSA_KEYGEN_H
