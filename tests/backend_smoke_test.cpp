#include "backend.h"
#include <cstdio>
using namespace emp;
using namespace mldsa;
#ifndef NP
#define NP 2
#endif
constexpr int nP = NP;
int main(int argc, char** argv) {
  int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));
  Backend<nP> bk(party, peer_port(), &pool);
  // Draw 8 authenticated random bits, open them (checked).
  auto r = bk.gmw().random(8);
  auto bits = bk.gmw().open(r, "smoke-open");
  bk.gmw().finish();
  if (party == 1) {
    std::printf("party1 opened bits:");
    for (auto b : bits) std::printf(" %d", (int)b);
    std::printf("\n comm=%lld bytes\n", (long long)bk.comm_bytes());
  }
  return 0;
}
