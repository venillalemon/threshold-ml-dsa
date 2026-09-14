#include "backend.h"
#include <emp-tool/circuits/frontend/frontend.h>
#include <cstdio>
using namespace emp;
using namespace mldsa;
#ifndef NP
#define NP 2
#endif
constexpr int nP = NP;
// 8-bit ripple adder: out = a + b (mod 256), plus carry-less. Exercises AND depth.
static const auto& adder() {
  using U8 = UInt_T<RecordCtx, 8>;
  static const auto c = frontend::compile<U8, U8>([](RecordCtx& ctx, U8 a, U8 b) {
    (void)ctx; return a + b;
  });
  return c;
}
int main(int argc, char** argv) {
  int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));
  Backend<nP> bk(party, peer_port(), &pool);
  const auto& prog = adder().program();
  // 16 input bits (two u8). Draw shares, open to learn plaintext inputs.
  auto in = bk.gmw().random(16);
  auto ibits = bk.gmw().open(in, "gmw-in-open");
  auto out = bk.gmw().evaluate(prog, in);
  auto obits = bk.gmw().open(out, "gmw-out-open");
  bk.gmw().finish();
  if (party == 1) {
    unsigned a=0,b=0,o=0;
    for (int k=0;k<8;k++){a|=ibits[k]<<k; b|=ibits[8+k]<<k; o|=obits[k]<<k;}
    std::printf("a=%u b=%u a+b mod256=%u  gmw=%u : %s\n", a,b,(a+b)&255,o,
                ((a+b)&255)==o?"OK":"FAIL");
  }
  return 0;
}
