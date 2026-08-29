// src/protocol/sign.h — Pi_MLDSA Sign, single-shot variant (T = 1).
//
// Paper steps -> code, same order, no repetition:
//   1a. F_PrepSign                              prepsign() (same session)
//   1b. mu = H(tr,m), c = H(mu,w1)              h256_stub + sample_in_ball_stub (CHEAT)
//   1c. <r0>_2 = <w0>_2 - c<e>_2,               Boolean circuit: c*e adder trees,
//       reveal [||r0||_inf < g2 - beta] only    one aggregated public bit
//       -> on reject return (c, bot, bot)
//   3.  z = y + c*s, ||z||_inf < g1 - beta      Boolean circuit (T=1: no index
//       -> on reject return (c, bot, bot)       selection), then open z publicly
//   4.  h = MakeHint(-c t0, Az - ct + c t0)     local on public values
//   5.  ||c t0||_inf >= g2 or HW(h) > omega     -> (c, z, bot), else (c, z, h)
//
// c*<e>_2 and c*<s>_2 exploit that c has TAU coefficients in {-1,+1}: each output
// coefficient is a signed sum of TAU (negacyclically indexed, possibly negated)
// small shares — ripple adder trees at CE_WIDTH bits, no general multiplier.
// Prefix sums stay in [-BETA, BETA] (m terms are bounded by m*ETA <= TAU*ETA).
#ifndef MLDSA_SIGN_H
#define MLDSA_SIGN_H

#include "circuit.h"
#include "keygen.h"
#include "prepsign.h"
#include "ref.h"
#include <emp-ag/emp-ag.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace mldsa {

constexpr int ce_bits(int bound) { // smallest W with |v| <= bound representable signed
  int w = 2;
  while ((int64_t(1) << (w - 1)) <= bound)
    ++w;
  return w;
}
constexpr int CE_WIDTH = ce_bits(BETA); // c*e / c*s accumulator width
constexpr int Z_WIDTH = 21;             // z = y + c*s, |z| <= GAMMA1 + BETA
static_assert(GAMMA1 + BETA < (1 << (Z_WIDTH - 1)), "Z_WIDTH too small");
static_assert(G2 + BETA < (1 << (OW0 - 1)), "w0 - c*e must fit OW0 bits signed");

// (c, z, h) with explicit bot markers: z_ok=false => (c,bot,bot) [r0 or z
// rejected]; h_ok=false => (c,z,bot). All public — this IS the signature.
struct Signature {
  std::vector<int32_t> c; // challenge polynomial, N coefficients in {-1,0,1}
  bool z_ok = false;
  std::vector<int32_t> z; // N*ELL coefficients, |z| <= GAMMA1-BETA (iff z_ok)
  bool h_ok = false;
  std::vector<uint8_t> h; // N*K hint bits (iff h_ok)
};

// c * <x>_2 for a length-`polys` vector of eta-small polynomials (x = e with
// polys = K, x = s with polys = ELL). Adopts the EW-bit two's-complement
// shares into the circuit and returns one CE_WIDTH-bit signed accumulator per
// coefficient. c is public and sparse with coefficients in {-1,+1}, so output
// coefficient i of poly p is a signed sum of |c_nz| negacyclically indexed,
// possibly negated shares — ripple adder trees, no general multiplier. Every
// negation is ~x plus the adder carry-in; prefix sums stay within
// [-BETA, BETA] (m terms are bounded by m*ETA), so CE_WIDTH never overflows.
template <int nP>
inline std::vector<std::array<emp::Bit_T<typename emp::AGMPCSession<nP>::ctx_t>, CE_WIDTH>>
cmul_shared(emp::AGMPCSession<nP>& sess, emp::ag::AShareBundleVec<nP>& x_2,
            const std::vector<std::pair<int, int>>& c_nz, int polys) {
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using Bit = emp::Bit_T<Ctx>;
  constexpr int EW = ETA == 2 ? 3 : 4;

  std::vector<std::array<Bit, CE_WIDTH>> x_bits((size_t)polys * N);
  for (int i = 0; i < polys * N; ++i) {
    auto xi =
        sess.template adopt_authenticated_input<emp::UInt_T<Ctx, EW>>(&x_2[(size_t)EW * i]);
    for (int k = 0; k < EW; ++k)
      x_bits[i][k] = xi[k];
    for (int k = EW; k < CE_WIDTH; ++k)
      x_bits[i][k] = xi[EW - 1]; // sign extend, free
  }

  std::vector<std::array<Bit, CE_WIDTH>> out((size_t)polys * N);
  for (int p = 0; p < polys; ++p)
    for (int i = 0; i < N; ++i) {
      Bit* acc = out[(size_t)p * N + i].data();
      bool first = true;
      for (const auto& [deg, sgn0] : c_nz) {
        int d = i - deg, sgn = sgn0;
        if (d < 0) { // negacyclic wrap: x^N = -1
          d += N;
          sgn = -sgn;
        }
        const Bit* src = x_bits[(size_t)p * N + d].data();
        if (first) {
          if (sgn > 0)
            for (int k = 0; k < CE_WIDTH; ++k)
              acc[k] = src[k];
          else { // -x = ~x + 1
            Bit inv[CE_WIDTH];
            for (int k = 0; k < CE_WIDTH; ++k)
              inv[k] = !src[k];
            add_const_ripple<CE_WIDTH>(sess.ctx(), acc, inv, 0, true);
          }
          first = false;
        } else {
          Bit term[CE_WIDTH];
          for (int k = 0; k < CE_WIDTH; ++k)
            term[k] = sgn < 0 ? !src[k] : src[k];
          add_ripple<CE_WIDTH>(sess.ctx(), acc, acc, term, sgn < 0);
        }
      }
    }
  return out;
}

// ---- Sign -------------------------------------------------------------------
template <int nP>
inline Signature sign(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer,
                      KeyPair<nP>& kp, const std::vector<uint32_t>& msg) {
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using Bit = emp::Bit_T<Ctx>;

  Signature sig;

#ifdef TEST
  auto lap_start = emp::clock_start();
  auto timer_lap = [&](const char* label) {
    if (party == 1)
      std::printf("  [timer] %-20s %8.1f ms\n", label, emp::time_from(lap_start) / 1000.0);
    lap_start = emp::clock_start();
  };
#endif

  // 1a. F_PrepSign -> (<w0>_2 wires, w1 public, <y>_2 shares), T = 1.
  std::vector<emp::UInt_T<Ctx, OW0>> w0_2;
  emp::ag::AShareBundleVec<nP> y_2;
  std::vector<uint32_t> w1;
  prepsign<nP>(sess, party, dealer, kp.A, w0_2, y_2, w1);
#ifdef TEST
  timer_lap("prepsign");
#endif

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

  // 1c. <r0>_2 = <w0>_2 - c*<e>_2; reveal only the STRICT predicate
  // [||r0||_inf < g2 - beta]. |w0 - c*e| <= g2 + beta < 2^(OW0-1), so
  // d := w0 - c*e is exact at OW0 bits, and the mod± 2g2 reduction never has
  // to be materialized: any d that would wrap (|d| > g2) reduces to
  // |r0| >= g2 - beta, which the strict bound rejects anyway. Hence
  // pass  <=>  -(g2-beta) < d < g2-beta  on the unreduced value.
  const auto ce = cmul_shared<nP>(sess, kp.e.two_share, c_nz, K);
  std::vector<Bit> r0_bad;
  r0_bad.reserve(COEFF_COUNT);
  for (int idx = 0; idx < COEFF_COUNT; ++idx) {
    Bit w0b[OW0];
    unpack_bits<OW0>(w0_2[idx], w0b);
    Bit nce[OW0]; // -ce sign-extended: ~ce, +1 via the adder carry-in
    for (int k = 0; k < OW0; ++k)
      nce[k] = !(k < CE_WIDTH ? ce[idx][k] : ce[idx][CE_WIDTH - 1]);
    Bit db[OW0];
    add_ripple<OW0>(sess.ctx(), db, w0b, nce, true); // d = w0 - ce
    const Bit in_lo = ge_const<OW0>(sess.ctx(), db, -(int64_t)(G2 - BETA) + 1); // d > -(g2-beta)
    const Bit in_hi = !ge_const<OW0>(sess.ctx(), db, (int64_t)(G2 - BETA));     // d < g2-beta
    r0_bad.push_back(!(in_lo & in_hi));
  }
  const Bit r0_reject = or_tree(r0_bad.data(), (int)r0_bad.size());
#ifdef TEST
  // Plaintext oracle: recompute every r0 bad bit from the opened e (keygen)
  // and w (prepsign) and compare with the circuit bit by bit.
  const auto centered = [](uint32_t x) {
    return x > (uint32_t)Q / 2 ? (int32_t)x - Q : (int32_t)x;
  };
  int r0_bad_count = 0, r0_mismatch = 0;
  for (int p = 0; p < K; ++p) {
    int64_t e_row[N], ce_row[N] = {0};
    for (int i = 0; i < N; ++i)
      e_row[i] = centered(kp.opened_e[p * N + i]);
    cpoly_mul_acc(sig.c, e_row, ce_row);
    for (int i = 0; i < N; ++i) {
      const int idx = p * N + i;
      int64_t d = ref_lowbits((int32_t)test_opened_w[idx], G2) - ce_row[i];
      if (d > G2) // mod± 2g2 into (-g2, g2]
        d -= 2 * (int64_t)G2;
      else if (d <= -G2)
        d += 2 * (int64_t)G2;
      const bool expect_bad = d >= G2 - BETA || d <= -(G2 - BETA); // strict acceptance
      const bool circuit_bad = (bool)sess.reveal(r0_bad[idx], emp::PUBLIC).value();
      r0_bad_count += circuit_bad;
      r0_mismatch += circuit_bad != expect_bad;
    }
  }
  if (party == 1)
    std::printf("  Sign: r0 out-of-range coefficients %d / %d, circuit/plaintext mismatches %d\n",
                r0_bad_count, COEFF_COUNT, r0_mismatch);
  emp::expecting(r0_mismatch == 0, "Sign test: r0 circuit disagrees with plaintext");
  timer_lap("r0 circuit");
#endif
  const auto r0_rej = sess.reveal(r0_reject, emp::PUBLIC);
  emp::expecting(r0_rej.has_value(), "Sign: r0 predicate reveal failed");
  if (r0_rej.value())
    return sig; // (c, bot, bot)

  // 3. z = y + c*s with the STRICT check ||z||_inf < g1 - beta; T = 1, so the
  // "first index j*" selection degenerates to this single check. Open z only
  // if it passes.
  const auto cs = cmul_shared<nP>(sess, kp.s.two_share, c_nz, ELL);
  std::vector<emp::UInt_T<Ctx, Z_WIDTH>> z_wires;
  z_wires.reserve(Y_COEFF_COUNT);
  std::vector<Bit> z_bad;
  z_bad.reserve(Y_COEFF_COUNT);
  for (int idx = 0; idx < Y_COEFF_COUNT; ++idx) {
    // <y>_2 stores y-1 (edabits.h affine encoding): z = (y-1) + cs + 1.
    auto yi = sess.template adopt_authenticated_input<emp::UInt_T<Ctx, Y_WIDTH>>(
        &y_2[(size_t)Y_WIDTH * idx]);
    Bit yb[Z_WIDTH], csb[Z_WIDTH], zb[Z_WIDTH];
    for (int k = 0; k < Z_WIDTH; ++k) {
      yb[k] = yi[k < Y_WIDTH ? k : Y_WIDTH - 1];
      csb[k] = cs[idx][k < CE_WIDTH ? k : CE_WIDTH - 1];
    }
    add_ripple<Z_WIDTH>(sess.ctx(), zb, yb, csb, true);
    const Bit ok = ge_const<Z_WIDTH>(sess.ctx(), zb, -(int64_t)(GAMMA1 - BETA) + 1) &
                   !ge_const<Z_WIDTH>(sess.ctx(), zb, (int64_t)(GAMMA1 - BETA));
    z_bad.push_back(!ok);
    z_wires.push_back(pack_bits<Z_WIDTH>(sess.ctx(), zb));
  }
  const Bit z_reject = or_tree(z_bad.data(), (int)z_bad.size());
#ifdef TEST
  // Plaintext oracle for z = y + c*s, bit-compared like r0 above.
  std::vector<int64_t> z_plain(Y_COEFF_COUNT);
  int z_bad_count = 0, z_mismatch = 0;
  for (int p = 0; p < ELL; ++p) {
    int64_t s_row[N], cs_row[N] = {0};
    for (int i = 0; i < N; ++i)
      s_row[i] = centered(kp.opened_s[p * N + i]);
    cpoly_mul_acc(sig.c, s_row, cs_row);
    for (int i = 0; i < N; ++i) {
      const int idx = p * N + i;
      z_plain[idx] = centered(test_opened_y[idx]) + cs_row[i];
      const bool expect_bad =
          z_plain[idx] >= GAMMA1 - BETA || z_plain[idx] <= -(GAMMA1 - BETA); // strict acceptance
      const bool circuit_bad = (bool)sess.reveal(z_bad[idx], emp::PUBLIC).value();
      z_bad_count += circuit_bad;
      z_mismatch += circuit_bad != expect_bad;
    }
  }
  if (party == 1)
    std::printf("  Sign: z out-of-range coefficients %d / %d, circuit/plaintext mismatches %d\n",
                z_bad_count, Y_COEFF_COUNT, z_mismatch);
  emp::expecting(z_mismatch == 0, "Sign test: z circuit disagrees with plaintext");
  timer_lap("z circuit");
#endif
  const auto z_rej = sess.reveal(z_reject, emp::PUBLIC);
  emp::expecting(z_rej.has_value(), "Sign: z predicate reveal failed");
  if (z_rej.value())
    return sig; // (c, bot, bot)

  sig.z.resize(Y_COEFF_COUNT);
  for (int i = 0; i < Y_COEFF_COUNT; ++i) {
    const auto zi = sess.reveal(z_wires[i], emp::PUBLIC);
    emp::expecting(zi.has_value(), "Sign: public z reveal failed");
    const int32_t raw = (int32_t)zi.value();
    sig.z[i] = raw - ((raw >> (Z_WIDTH - 1)) << Z_WIDTH); // sign extend
#ifdef TEST
    emp::expecting(sig.z[i] == (int32_t)z_plain[i], "Sign test: opened z != plaintext y + c*s");
#endif
  }
  sig.z_ok = true;
#ifdef TEST
  timer_lap("z reveal");
#endif

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
      const uint32_t u = modq_i64(az - ct[i] + ct0[i]); // Az - ct + c*t0
      const uint8_t hi = (uint8_t)ref_makehint(-ct0[i], (int32_t)u, G2);
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
