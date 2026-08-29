// main.cpp — run Pi_MLDSA KeyGen, then a single-shot (T=1) Sign.
#include "keygen.h"
#include "sign.h"
#include <emp-ag/emp-ag.h>

#include <cstdint>
#include <cstdio>
#include <random>

using namespace emp;
using namespace mldsa;

#ifndef MLDSA_NP
#define MLDSA_NP 3
#endif
constexpr int nP = MLDSA_NP;

int main(int argc, char** argv) {
  const int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));

  const auto setup_start = clock_start();
  AGMPCSession<nP> sess(party, peer_port(), /*ip table*/ nullptr, &pool);
  const double setup_ms = time_from(setup_start) / 1000.0;

  // Demo-only SPDZ offline from a common seed (spdz.h FakeDealer): the full
  // MAC key lives only inside the dealer; protocols see my_alpha and deal().
  FakeDealer<nP> dealer(party, 777);

  const auto keygen_start = clock_start();
  KeyPair<nP> kp = keygen<nP>(sess, party, dealer);
  const double keygen_ms = time_from(keygen_start) / 1000.0;

  expecting((int)kp.t1.size() == COEFF_COUNT, "KeyGen: wrong t1 length");
  expecting((int)kp.s.q_share.size() == Y_COEFF_COUNT, "KeyGen: wrong s length");
  expecting((int)kp.e.q_share.size() == COEFF_COUNT, "KeyGen: wrong e length");

  const auto sign_start = clock_start();
  const std::vector<uint32_t> msg = {0x48656c6cu, 0x6f204d4cu, 0x2d445341u}; // "Hello ML-DSA"
  Signature sig = mldsa::sign<nP>(sess, party, dealer, kp, msg);
  const double sign_ms = time_from(sign_start) / 1000.0;

  expecting((int)sig.c.size() == N, "Sign: wrong c length");
  if (sig.z_ok)
    expecting((int)sig.z.size() == Y_COEFF_COUNT, "Sign: wrong z length");
  if (sig.h_ok)
    expecting((int)sig.h.size() == COEFF_COUNT, "Sign: wrong h length");

  if (party == 1) {
    std::printf("  pk: rho[0]=%08x, t1=%zu coefficients | tr[0]=%08x\n", kp.rho[0], kp.t1.size(),
                kp.tr[0]);
    if (!sig.z_ok) {
      std::printf(
          "\033[31m  signature: (c, bot, bot) — rejected, retry with fresh PrepSign (T=1)\033[0m\n");
    } else if (!sig.h_ok) {
      std::printf("\033[31m  signature: (c, z, bot) — hint checks failed\033[0m\n");
    } else {
      int hw = 0, z_max = 0;
      for (uint8_t b : sig.h)
        hw += b;
      for (int32_t v : sig.z)
        z_max = std::abs(v) > z_max ? std::abs(v) : z_max;
      std::printf("\033[32m  signature: (c, z, h) — ||z||_inf=%d (< %d), HW(h)=%d (<= %d)\033[0m\n",
                  z_max, GAMMA1 - BETA, hw, OMEGA);
    }
    std::printf("  setup %.1f ms | KeyGen %.1f ms | Sign %.1f ms | ANDs garbled %llu\n", setup_ms,
                keygen_ms, sign_ms, (unsigned long long)sess.num_and());
    std::printf("  KeyGen + Sign complete\n");
  }
  return 0;
}
