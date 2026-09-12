// main.cpp — run Pi_MLDSA KeyGen, then a single-shot (T=1) Sign.
#include "keygen.h"
// MLDSA_SLOT selects the post-challenge design: 1 = slot-restricted boundary
// test (threshold-mldsa-slot.pdf), 0 = the two-round per-coefficient baseline.
#if MLDSA_SLOT
#include "sign.h"
#else
#include "sign_2round.h"
#endif
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

// This party's total bytes on the wire (sent + received, all meshes).
static int64_t comm_bytes(AGMPCSession<nP>& sess) {
  if (auto* grp = sess.group()) {
    int64_t total = 0;
    for (int m = 0; m < grp->size(); ++m)
      total += grp->mesh(m).count();
    return total;
  }
  return sess.io().count();
}

// NOTE on "communication rounds": NetIOMP gives every peer pair TWO
// unidirectional sockets (ios[p] send-only, ios2[p] recv-only -- see
// backend/netmp.h send_data/recv_data), so IOChannel::rounds (a same-channel
// send/recv direction-switch counter) can never move past the initial
// handshake on this transport and is not a usable round-complexity proxy
// here. Instead, sign.h counts actual reveal()/decode barriers -- see
// `online_reveal_rounds` below -- which is what ChunkedSession::reveal
// treats as a synchronization point (flush pending gates, then decode).

int main(int argc, char** argv) {
  const int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));

  const auto setup_start = clock_start();
  AGMPCSession<nP> sess(party, peer_port(), /*ip table*/ nullptr, &pool);
  const double setup_ms = time_from(setup_start) / 1000.0;
  const int64_t setup_comm = comm_bytes(sess);

  // Demo-only offline from a common seed (spdz.h FakeDealer): every party's
  // BDOZ slope lives only inside the dealer; protocols see their own slope.
  FakeDealer<nP> dealer(party, 777);

  const auto keygen_start = clock_start();
  KeyPair<nP> kp = keygen<nP>(sess, party, dealer);
  const double keygen_ms = time_from(keygen_start) / 1000.0;
  const uint64_t keygen_ands = sess.num_and();
  const uint64_t keygen_edabit_ands = g_edabit_ands; // exclude correlation gates from G1
  const int64_t keygen_comm = comm_bytes(sess) - setup_comm;

  expecting((int)kp.t1.size() == COEFF_COUNT, "KeyGen: wrong t1 length");
  expecting((int)kp.s.fq_share.size() == Y_COEFF_COUNT, "KeyGen: wrong s length");
  expecting((int)kp.e.fq_share.size() == COEFF_COUNT, "KeyGen: wrong e length");

  // Sign = F_PrepSign (offline, independent of msg) + the online tail (needs
  // msg to derive the challenge c). sign()'s after_prepsign hook fires right
  // at that boundary so we can split Sign's cost the same way KeyGen/Sign are
  // split above -- everything from here to the hook is "offline-able" prep;
  // everything after is genuinely online (latency the signer pays once the
  // message shows up). online_reveal_rounds counts that tail's actual
  // reveal() barriers (see sign.h) -- the real round-complexity number.
  const int64_t pre_sign_comm = comm_bytes(sess);
  double prepsign_ms = 0.0;
  uint64_t prepsign_ands = 0;
  int64_t prepsign_comm = 0;
  int64_t online_reveal_rounds = 0, g2_ands = 0;
  const auto sign_start = clock_start();
  const std::vector<uint32_t> msg = {0x48656c6cu, 0x6f204d4cu, 0x2d445341u}; // "Hello ML-DSA"
  Signature sig = mldsa::sign<nP>(
      sess, party, dealer, kp, msg,
      [&] {
        prepsign_ms = time_from(sign_start) / 1000.0;
        prepsign_ands = sess.num_and() - keygen_ands;
        prepsign_comm = comm_bytes(sess) - pre_sign_comm;
      },
      &online_reveal_rounds, &g2_ands);
  const double sign_ms = time_from(sign_start) / 1000.0;
  const uint64_t sign_ands = sess.num_and() - keygen_ands;
  const int64_t sign_comm = comm_bytes(sess) - setup_comm - keygen_comm;
  // Online tail = everything sign() did after the hook fired.
  const double online_ms = sign_ms - prepsign_ms;
  const uint64_t online_ands = sign_ands - prepsign_ands;
  const int64_t online_comm = sign_comm - prepsign_comm;

  expecting((int)sig.c.size() == N, "Sign: wrong c length");
  if (sig.z_ok)
    expecting((int)sig.z.size() == Y_COEFF_COUNT, "Sign: wrong z length");
  if (sig.h_ok)
    expecting((int)sig.h.size() == COEFF_COUNT, "Sign: wrong h length");

  // Every party: its own online sent+received bytes (network total = sum / 2).
  const int64_t sign_edabit_ands = (int64_t)(g_edabit_ands - keygen_edabit_ands);
  std::printf("COMM v=%d party=%d online_comm=%lld prepsign_comm=%lld g1_ands=%lld g2_ands=%lld "
              "edabit_ands=%lld\n",
              MLDSA_SLOT, party, (long long)online_comm, (long long)prepsign_comm,
              (long long)((int64_t)sign_ands - g2_ands - sign_edabit_ands), (long long)g2_ands,
              (long long)sign_edabit_ands);
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
    std::printf("  KeyGen %7.1f ms | comm %6.1f MB | ANDs %llu\n", keygen_ms, keygen_comm / 1e6,
                (unsigned long long)keygen_ands);
    std::printf("  Sign   %7.1f ms | comm %6.1f MB | ANDs %llu\n", sign_ms, sign_comm / 1e6,
                (unsigned long long)sign_ands);
    std::printf("    - PrepSign (offline, no msg) %7.1f ms | comm %6.1f MB | ANDs %llu\n", prepsign_ms,
                prepsign_comm / 1e6, (unsigned long long)prepsign_ands);
    std::printf("    - online   (after msg known) %7.1f ms | comm %6.1f MB | ANDs %llu | reveal "
                "rounds %lld\n",
                online_ms, online_comm / 1e6, (unsigned long long)online_ands,
                (long long)online_reveal_rounds);
    // Machine-readable line for benchmark scripts (bench.sh greps "STATS").
    std::printf("STATS np=%d param=%s setup_ms=%.1f keygen_ms=%.1f sign_ms=%.1f "
                "prepsign_ms=%.1f online_ms=%.1f "
                "keygen_ands=%llu sign_ands=%llu prepsign_ands=%llu online_ands=%llu "
                "setup_comm=%lld keygen_comm=%lld sign_comm=%lld prepsign_comm=%lld online_comm=%lld "
                "online_reveal_rounds=%lld accepted=%d\n",
                nP, param_name(PARAM), setup_ms, keygen_ms, sign_ms, prepsign_ms, online_ms,
                (unsigned long long)keygen_ands, (unsigned long long)sign_ands,
                (unsigned long long)prepsign_ands, (unsigned long long)online_ands,
                (long long)setup_comm, (long long)keygen_comm, (long long)sign_comm,
                (long long)prepsign_comm, (long long)online_comm, (long long)online_reveal_rounds,
                sig.z_ok && sig.h_ok ? 1 : 0);
    std::printf("  KeyGen + Sign complete\n");
  }
  return 0;
}
