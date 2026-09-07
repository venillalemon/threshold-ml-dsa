// src/protocol/sign.h — Pi_MLDSA Sign, single-shot variant (T = 1), online-minimal.
//
// Paper steps -> code:
//   1a. F_PrepSign                              prepsign()            [offline]
//   1a'. edaBit masks for z and r0; garble      Fq_edabits + sess.prepare()
//        the online circuit                                            [offline]
//   1b. mu = H(tr,m), c = H(mu,w1)              h256_stub + sample_in_ball_stub (CHEAT)
//   1c+3. <z>_q = <y>_q + c<s>_q,               LOCAL (linear on SPDZ shares)
//         <d>_q = <w>_q - w1*2g2 - c<e>_q       LOCAL (= w0 - c*e, see below)
//         open z + r_z, d + r_d                 ONE all-to-all open + MACCheck
//         circuit: z = (cz - r_z) mod q,        ONE prepared circuit run:
//         d = (cd - r_d) mod q, both norm        output = (accept, accept ? z : 0)
//         predicates -> accept, z via MUX        ONE batched reveal
//       -> on reject return (c, bot, bot)
//   4.  h = MakeHint(-c t0, Az - ct + c t0)     local on public values
//   5.  ||c t0||_inf >= g2 or HW(h) > omega     -> (c, z, bot), else (c, z, h)
//
// c never touches a Boolean circuit: c*s and c*e are computed on the SPDZ
// arithmetic shares (free), and w0 - c*e is obtained as <w>_q - w1*2*gamma2 -
// c<e>_q, which equals FIPS Decompose's w0 minus c*e mod q (the identity
// w = w1*2g2 + w0 holds mod q including the top-of-range corner). Only the
// two range predicates are non-linear; they run on a CHALLENGE-INDEPENDENT
// prepared circuit whose inputs are the edaBit masks (soldered) and the two
// masked openings (public). Online rounds: open (1) + MACCheck (1) +
// run(prepared) (~4) + reveal (1). Both predicates fold into one accept bit
// and z is released through MUX(accept, z, 0), so a rejected run reveals
// nothing but the bit.
#ifndef MLDSA_SIGN_H
#define MLDSA_SIGN_H

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

constexpr int ZW = Y_COEFF_COUNT * L;  // z / r_z / cz, L-bit mod-q representatives
constexpr int DW = COEFF_COUNT * L;    // d / r_d / cd
constexpr int OUTW = 1 + ZW;           // accept || z (mod-q form, gated)
static_assert(GAMMA1 + BETA < Q / 2 && G2 + BETA < Q / 2, "centered bounds must be < q/2");

template <class Ctx> using ZV = emp::BitVec_T<Ctx, ZW>;
template <class Ctx> using DV = emp::BitVec_T<Ctx, DW>;
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

// ---- local SPDZ arithmetic --------------------------------------------------

// acc[N] += c * x over R_q on authenticated shares; c sparse in {-1,0,1}
// given as (degree, sign) pairs. Purely local (linear).
inline void cpoly_mul_auth(const std::vector<std::pair<int, int>>& c_nz, const AuthShare* x,
                           AuthShare* acc) {
  for (const auto& [deg, sgn0] : c_nz)
    for (int i = 0; i < N; ++i) {
      int d = i - deg, sgn = sgn0;
      if (d < 0) { // negacyclic wrap: x^N = -1
        d += N;
        sgn = -sgn;
      }
      acc[i] = acc[i] + (sgn > 0 ? x[d] : x[d] * (uint32_t)(Q - 1));
    }
}

// ---- the challenge-independent online circuit ------------------------------

// v = (c - r) mod q with c public, r shared, both in [0, q). Two L-bit ripple
// chains (2L AND): d = c - r + 2^L, then +Q iff borrowed. Q's bits are
// compile-time constants, so the conditional +Q chain needs no extra AND.
template <class Ctx>
inline void sub_modq_bits(Ctx& ctx, Bit_T<Ctx>* v, const Bit_T<Ctx>* c, const Bit_T<Ctx>* r) {
  using Bit = Bit_T<Ctx>;
  Bit d[L], carry = Bit::constant(ctx, true);
  for (int i = 0; i < L; ++i)
    fa(d[i], carry, c[i], !r[i], carry);
  const Bit need = !carry; // c < r
  carry = Bit::constant(ctx, false);
  for (int i = 0; i < L; ++i) {
    if (((Q >> i) & 1) != 0)
      fa(v[i], carry, d[i], need, carry);
    else
      ha(v[i], carry, d[i], carry);
  }
}

// Strict centered bound |v_centered| < B for v in [0, q): v < B or v > q - B.
template <class Ctx>
inline Bit_T<Ctx> in_centered_bound(Ctx& ctx, const Bit_T<Ctx>* v, int64_t B) {
  return !ge_const_u<L>(ctx, v, (uint64_t)B) | ge_const_u<L>(ctx, v, (uint64_t)(Q - B + 1));
}

// Output bit 0 is accept; bits 1.. are z (mod-q representatives) gated by accept.
template <class Ctx>
inline OUTV<Ctx> online_body(Ctx& ctx, const ZV<Ctx>& rz, const DV<Ctx>& rd, const ZV<Ctx>& cz,
                             const DV<Ctx>& cd) {
  using Bit = Bit_T<Ctx>;
  std::vector<Bit> bad;
  bad.reserve((size_t)COEFF_COUNT + Y_COEFF_COUNT);
  std::vector<Bit> zb((size_t)ZW);
  Bit c[L], r[L], v[L];

  for (int idx = 0; idx < Y_COEFF_COUNT; ++idx) {
    for (int k = 0; k < L; ++k) {
      c[k] = cz[idx * L + k];
      r[k] = rz[idx * L + k];
    }
    sub_modq_bits(ctx, v, c, r);
    bad.push_back(!in_centered_bound(ctx, v, (int64_t)(GAMMA1 - BETA)));
    for (int k = 0; k < L; ++k)
      zb[(size_t)(idx * L + k)] = v[k];
  }
  for (int idx = 0; idx < COEFF_COUNT; ++idx) {
    for (int k = 0; k < L; ++k) {
      c[k] = cd[idx * L + k];
      r[k] = rd[idx * L + k];
    }
    sub_modq_bits(ctx, v, c, r);
    bad.push_back(!in_centered_bound(ctx, v, (int64_t)(G2 - BETA)));
  }

  const Bit accept = !or_tree(bad.data(), (int)bad.size());
  std::vector<Bit> outb((size_t)OUTW);
  outb[0] = accept;
  for (int b = 0; b < ZW; ++b)
    outb[(size_t)(1 + b)] = accept & zb[(size_t)b]; // MUX(accept, z, 0)
  return OUTV<Ctx>::from_bit_values(ctx, outb.data());
}

// Compiled once per process: the circuit depends on the parameter set only.
inline const auto& online_circuit() {
  using emp::RecordCtx;
  static const auto c = emp::frontend::compile<ZV<RecordCtx>, DV<RecordCtx>, ZV<RecordCtx>,
                                               DV<RecordCtx>>(
      [](RecordCtx& ctx, ZV<RecordCtx> rz, DV<RecordCtx> rd, ZV<RecordCtx> cz, DV<RecordCtx> cd) {
        return online_body(ctx, rz, rd, cz, cd);
      });
  return c;
}

// No-op default for sign()'s after_prepsign hook (see below).
struct NoopHook {
  void operator()() const {}
};

// ---- Sign -------------------------------------------------------------------
// `after_prepsign` fires once, right after the offline part (F_PrepSign, the
// edaBit masks and the prepared garbling) and before any message-dependent
// step -- lets a caller snapshot comm counters at exactly the offline/online
// boundary. `online_rounds_out`, if non-null, receives the number of online
// synchronization barriers (open+MACCheck, run(prepared), reveal = 3).
template <int nP, class Hook = NoopHook>
inline Signature sign(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer,
                      KeyPair<nP>& kp, const std::vector<uint32_t>& msg,
                      Hook after_prepsign = Hook{}, int64_t* online_rounds_out = nullptr) {
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using Wire = typename Ctx::Wire;
  using UL = emp::UInt_T<Ctx, L>;

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
  // 1a. F_PrepSign -> (<y>_q, <w>_q, w1 public), T = 1.
  std::vector<AuthShare> y_q, w_q;
  std::vector<uint32_t> w1;
  prepsign<nP>(sess, party, dealer, kp.A, y_q, w_q, w1);
#ifdef TEST
  timer_lap("prepsign");
#endif

  // 1a'. edaBit masks for the two A2B conversions, then garble the
  // (challenge-independent) online circuit: OTs, triples, rows and the COT
  // check all finish before the message exists.
  SharePair<nP, L, true> rz = Fq_edabits<nP, L>(sess, party, dealer, Y_COEFF_COUNT);
  SharePair<nP, L, true> rd = Fq_edabits<nP, L>(sess, party, dealer, COEFF_COUNT);
  auto prepared = sess.prepare(online_circuit());
#ifdef TEST
  timer_lap("edabits + prepare");
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

  // 1c+3 (linear part, local): <z>_q = <y>_q + c<s>_q,
  // <d>_q = <w>_q - w1*2g2 - c<e>_q  (= w0 - c*e mod q).
  std::vector<AuthShare> z_q = y_q;
  for (int p = 0; p < ELL; ++p)
    cpoly_mul_auth(c_nz, &kp.s.q_share[(size_t)p * N], &z_q[(size_t)p * N]);
  std::vector<AuthShare> d_q = w_q;
  for (int idx = 0; idx < COEFF_COUNT; ++idx)
    d_q[(size_t)idx] = add_public(d_q[(size_t)idx],
                                  modq_i64(-(int64_t)w1[(size_t)idx] * 2 * G2), party,
                                  dealer.my_alpha);
  {
    std::vector<AuthShare> ce((size_t)COEFF_COUNT);
    for (int p = 0; p < K; ++p)
      cpoly_mul_auth(c_nz, &kp.e.q_share[(size_t)p * N], &ce[(size_t)p * N]);
    for (int idx = 0; idx < COEFF_COUNT; ++idx)
      d_q[(size_t)idx] = d_q[(size_t)idx] + ce[(size_t)idx] * (uint32_t)(Q - 1);
  }

  // Mask-and-open (F_A2B step 2): cz = z + r_z, cd = d + r_d, one all-to-all
  // open for both vectors, one batched MACCheck.
  std::vector<uint32_t> open_val((size_t)Y_COEFF_COUNT + COEFF_COUNT), open_mac(open_val.size());
  for (int i = 0; i < Y_COEFF_COUNT; ++i) {
    const AuthShare m = z_q[(size_t)i] + rz.q_share[(size_t)i];
    open_val[(size_t)i] = m.val;
    open_mac[(size_t)i] = m.mac;
  }
  for (int i = 0; i < COEFF_COUNT; ++i) {
    const AuthShare m = d_q[(size_t)i] + rd.q_share[(size_t)i];
    open_val[(size_t)Y_COEFF_COUNT + i] = m.val;
    open_mac[(size_t)Y_COEFF_COUNT + i] = m.mac;
  }
#ifdef TAMPER_C
  if (party == 1)
    open_val[0] ^= 1;
#endif
  const std::vector<uint32_t> opened_masked = open_additive_modq(sess.io(), party, open_val);
  spdz_maccheck(sess.io(), party, dealer.my_alpha, open_mac, opened_masked);

  // Circuit inputs: the masks (soldered) and the openings (public).
  std::vector<Wire> rzw((size_t)ZW), rdw((size_t)DW);
  for (int i = 0; i < Y_COEFF_COUNT; ++i)
    sess.template adopt_authenticated_input<UL>(&rz.two_share[(size_t)L * i])
        .pack_wires(&rzw[(size_t)i * L]);
  for (int i = 0; i < COEFF_COUNT; ++i)
    sess.template adopt_authenticated_input<UL>(&rd.two_share[(size_t)L * i])
        .pack_wires(&rdw[(size_t)i * L]);
  const ZV<Ctx> rz_v = ZV<Ctx>::from_wires(sess.ctx(), rzw.data());
  const DV<Ctx> rd_v = DV<Ctx>::from_wires(sess.ctx(), rdw.data());
  typename ZV<Ctx>::clear_t czb{};
  typename DV<Ctx>::clear_t cdb{};
  for (int i = 0; i < Y_COEFF_COUNT; ++i)
    for (int k = 0; k < L; ++k)
      czb[(size_t)(i * L + k)] = ((opened_masked[(size_t)i] >> k) & 1) != 0;
  for (int i = 0; i < COEFF_COUNT; ++i)
    for (int k = 0; k < L; ++k)
      cdb[(size_t)(i * L + k)] = ((opened_masked[(size_t)Y_COEFF_COUNT + i] >> k) & 1) != 0;
  const ZV<Ctx> cz_v = sess.template input<ZV<Ctx>>(emp::PUBLIC, czb);
  const DV<Ctx> cd_v = sess.template input<DV<Ctx>>(emp::PUBLIC, cdb);

  // One prepared run, one batched reveal.
  const OUTV<Ctx> out = sess.run(std::move(prepared), rz_v, rd_v, cz_v, cd_v);
  const auto opened = sess.reveal(out, emp::PUBLIC);
  emp::expecting(opened.has_value(), "Sign: online reveal failed");
  if (online_rounds_out)
    *online_rounds_out = 3;
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
  if (!accept) {
    for (int b = 1; b < OUTW; ++b)
      emp::expecting(!ob[(size_t)b], "Sign test: rejected run leaked a nonzero z bit");
  }
#endif

  if (!accept)
    return sig; // (c, bot, bot)

  sig.z.resize(Y_COEFF_COUNT);
  for (int i = 0; i < Y_COEFF_COUNT; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < L; ++k)
      v |= (uint32_t)ob[(size_t)(1 + i * L + k)] << k;
    sig.z[i] = v > (uint32_t)Q / 2 ? (int32_t)v - Q : (int32_t)v; // centre
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
