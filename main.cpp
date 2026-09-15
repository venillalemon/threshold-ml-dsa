// main.cpp — run Pi_MLDSA KeyGen, then a single-shot (T=1) Sign, on the
// mixed-security backend: authenticated GMW for the offline circuits and WRK17b
// four-row garbling (garbled offline) for the online circuit.
#include "backend.h"
#include "keygen.h"
#if MLDSA_SLOT
#include "sign_slot.h"
#else
#include "sign_2round.h"
#endif

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace emp;
using namespace mldsa;

#ifndef MLDSA_NP
#define MLDSA_NP 3
#endif
constexpr int nP = MLDSA_NP;

int main(int argc, char** argv) {
  const int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));

  // Demo-only offline (dealer.h FakeDealer): every correlation and every
  // authentication key is dealt from a common seed; protocols only ever see
  // their own slot. Seed from the shared port so all parties agree yet each
  // run differs; EMP_SEED reproduces a specific run.
  const char* seed_env = std::getenv("EMP_SEED");
  const uint32_t seed = seed_env ? (uint32_t)std::strtoul(seed_env, nullptr, 10)
                                 : (0x9e3779b9u ^ (uint32_t)peer_port());
  FakeDealer<nP> dealer(party, seed);

  // The backend takes the dealer's Delta so dealer-fabricated Boolean shares
  // authenticate under it (drop the argument for a private random key).
  const auto setup_start = clock_start();
  Backend<nP> bk(party, peer_port(), &pool, dealer.my_delta);
  const double setup_ms = time_from(setup_start) / 1000.0;
  const int64_t setup_comm = bk.comm_bytes();

  const auto keygen_start = clock_start();
  KeyPair<nP> kp = keygen<nP>(bk, party, dealer);
  const double keygen_ms = time_from(keygen_start) / 1000.0;
  const uint64_t keygen_ands = bk.gmw().ands_evaluated;
  const int64_t keygen_comm = bk.comm_bytes() - setup_comm;

  expecting((int)kp.t1.size() == COEFF_COUNT, "KeyGen: wrong t1 length");
  expecting((int)kp.s.fq_share.size() == Y_COEFF_COUNT, "KeyGen: wrong s length");
  expecting((int)kp.e.fq_share.size() == COEFF_COUNT, "KeyGen: wrong e length");

  // Sign = offline (edaBits, producers via GMW, C_post garbled + delivered) +
  // the online tail (one delta_H open + one flight of labels + local decode).
  // after_prepsign fires at the offline/online boundary.
  const int64_t pre_sign_comm = bk.comm_bytes();
  double offline_ms = 0.0;
  int64_t offline_comm = 0;
  int64_t online_reveal_rounds = 0, g2_ands = 0;
  const auto sign_start = clock_start();
  const std::vector<uint32_t> msg = {0x48656c6cu, 0x6f204d4cu, 0x2d445341u}; // "Hello ML-DSA"
  Signature sig = mldsa::sign<nP>(
      bk, party, dealer, kp, msg,
      [&] {
        offline_ms = time_from(sign_start) / 1000.0;
        offline_comm = bk.comm_bytes() - pre_sign_comm;
      },
      &online_reveal_rounds, &g2_ands);
  const double sign_ms = time_from(sign_start) / 1000.0;
  const int64_t sign_comm = bk.comm_bytes() - setup_comm - keygen_comm;
  const double online_ms = sign_ms - offline_ms;
  const int64_t online_comm = sign_comm - offline_comm;
  const uint64_t sign_ands = bk.gmw().ands_evaluated - keygen_ands;

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
    std::printf("  setup  %7.1f ms | comm %6.1f MB\n", setup_ms, setup_comm / 1e6);
    std::printf("  KeyGen %7.1f ms | comm %6.1f MB | GMW ANDs %llu\n", keygen_ms,
                keygen_comm / 1e6, (unsigned long long)keygen_ands);
    std::printf("  Sign   %7.1f ms | comm %6.1f MB | GMW ANDs %llu | G2(WRK) ANDs %lld\n", sign_ms,
                sign_comm / 1e6, (unsigned long long)sign_ands, (long long)g2_ands);
    std::printf("    - offline (no msg)   %7.1f ms | comm %6.1f MB\n", offline_ms,
                offline_comm / 1e6);
    std::printf("    - online  (after msg) %7.1f ms | comm %6.1f MB | flights %lld\n", online_ms,
                online_comm / 1e6, (long long)online_reveal_rounds);
    std::printf("STATS np=%d param=%s slot=%d setup_ms=%.1f keygen_ms=%.1f sign_ms=%.1f "
                "offline_ms=%.1f online_ms=%.1f keygen_ands=%llu sign_ands=%llu g2_ands=%lld "
                "setup_comm=%lld keygen_comm=%lld offline_comm=%lld online_comm=%lld "
                "online_flights=%lld accepted=%d\n",
                nP, param_name(PARAM), (int)MLDSA_SLOT, setup_ms, keygen_ms, sign_ms, offline_ms,
                online_ms, (unsigned long long)keygen_ands, (unsigned long long)sign_ands,
                (long long)g2_ands, (long long)setup_comm, (long long)keygen_comm,
                (long long)offline_comm, (long long)online_comm, (long long)online_reveal_rounds,
                sig.z_ok && sig.h_ok ? 1 : 0);
    std::printf("  KeyGen + Sign complete\n");
  }
  return 0;
}
