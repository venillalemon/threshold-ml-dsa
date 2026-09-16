// main.cpp — run Pi_MLDSA KeyGen, then a single-shot (T=1) Sign, on the
// mixed-security backend, and report gate counts and communication.
//
// Report (party 1 prints):
//   offline gates : ANDs of C_pre (the offline GMW circuit)
//   online  gates : ANDs of C_post (the WRK circuit evaluated after the message)
//   offline comm  : setup (COT) + KeyGen + everything before the message
//   online  comm  : everything after the message (flight 1 + flight 2)
// each communication figure in two accountings:
//   total = sum over ALL parties of bytes SENT (every byte counted once)
//   P1    = bytes sent + received by party 1 (the evaluator's own link)
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

  const char* seed_env = std::getenv("EMP_SEED");
  const uint32_t seed = seed_env ? (uint32_t)std::strtoul(seed_env, nullptr, 10)
                                 : (0x9e3779b9u ^ (uint32_t)peer_port());
  FakeDealer<nP> dealer(party, seed);
  Backend<nP> bk(party, peer_port(), &pool, dealer.my_delta);

  KeyPair<nP> kp = keygen<nP>(bk, party, dealer);
  expecting((int)kp.t1.size() == COEFF_COUNT, "KeyGen: wrong t1 length");

  int64_t off_sent = 0, off_recv = 0;
  int64_t rounds = 0, g2 = 0, g1 = 0;
  const std::vector<uint32_t> msg = {0x48656c6cu, 0x6f204d4cu, 0x2d445341u}; // "Hello ML-DSA"
  Signature sig = mldsa::sign<nP>(
      bk, party, dealer, kp, msg,
      [&] { // offline/online boundary: freeze the offline counters
        off_sent = bk.bytes_sent();
        off_recv = bk.bytes_recv();
      },
      &rounds, &g2, &g1);
  const int64_t on_sent = bk.bytes_sent() - off_sent;
  const int64_t on_recv = bk.bytes_recv() - off_recv;
  expecting((int)sig.c.size() == N, "Sign: wrong c length");

  // Gather every party's SENT bytes at P1 (after all measurements are frozen;
  // this exchange is not part of the protocol and is not counted).
  int64_t tot_off = off_sent, tot_on = on_sent;
  if (party != 1) {
    int64_t mine[2] = {off_sent, on_sent};
    bk.io().send_data(1, mine, sizeof(mine));
    bk.io().flush(1);
    return 0;
  }
  for (int p = 2; p <= nP; ++p) {
    int64_t theirs[2];
    bk.io().recv_data(p, theirs, sizeof(theirs));
    tot_off += theirs[0];
    tot_on += theirs[1];
  }

  std::printf("np=%d param=%s slot=%d online_flights=%lld accepted=%d\n", nP, param_name(PARAM),
              (int)MLDSA_SLOT, (long long)rounds, sig.z_ok && sig.h_ok ? 1 : 0);
  std::printf("offline gates (C_pre, GMW)  %12lld\n", (long long)g1);
  std::printf("online  gates (C_post, WRK) %12lld\n", (long long)g2);
  std::printf("offline comm  total          %12.3f MB\n", tot_off / 1e6);
  std::printf("offline comm  P1 sent+recv   %12.3f MB\n", (off_sent + off_recv) / 1e6);
  std::printf("online  comm  total          %12.3f MB\n", tot_on / 1e6);
  std::printf("online  comm  P1 sent+recv   %12.3f MB\n", (on_sent + on_recv) / 1e6);
  return 0;
}
