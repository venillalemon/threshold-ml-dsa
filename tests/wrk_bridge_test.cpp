#include "backend.h"
#include "wrk_phase2.h"
#include <emp-tool/circuits/frontend/frontend.h>
#include <cstdio>
using namespace emp;
using namespace mldsa;
#ifndef NP
#define NP 2
#endif
constexpr int nP = NP;

// C_post-like toy: input 0 = fixed secret x, input 1 = late public y.
// out0 = x AND y, out1 = x XOR y.  Exercises a product, a fixed and a late port.
static const auto& toy_circuit() {
  using B1 = BitVec_T<RecordCtx, 1>;
  using B2 = BitVec_T<RecordCtx, 2>;
  static const auto c = frontend::compile<B1, B1>([](RecordCtx& ctx, B1 x, B1 y) {
    Bit_T<RecordCtx> xb(ctx, 0), yb(ctx, 0);
    { typename RecordCtx::Wire w; x.pack_wires(&w); xb = Bit_T<RecordCtx>(ctx, w); }
    { typename RecordCtx::Wire w; y.pack_wires(&w); yb = Bit_T<RecordCtx>(ctx, w); }
    Bit_T<RecordCtx> outb[2] = { xb & yb, xb ^ yb };
    return B2::from_bit_values(ctx, outb);
  });
  return c;
}

int main(int argc, char** argv) {
  int party = parse_party(argv, nP);
  ThreadPool pool(ag::default_pool_size(nP));
  Backend<nP> bk(party, peer_port(), &pool);
  const auto& circ = toy_circuit();
  const auto& prog = circ.program();

  // Fixed secret x: draw one GMW share, open it so the test knows the truth.
  auto xs = bk.gmw().random(1);
  auto xbits = bk.gmw().open(xs, "toy-x-open-for-test");
  const uint8_t x = xbits[0];

  auto off = wrk_offline<nP>(bk, prog, xs, /*fixed_inputs=*/1);

  // Online with public y = 1.
  const uint8_t y = 1;
  auto out = wrk_online<nP>(bk, off, std::vector<uint8_t>{y});
  bk.gmw().finish();
  if (party == 1) {
    std::printf("x=%d y=%d -> out0=%d (want %d) out1=%d (want %d) : %s\n",
                x, y, out->at(0), x & y, out->at(1), x ^ y,
                (out->at(0) == (x & y) && out->at(1) == (x ^ y)) ? "OK" : "FAIL");
  }
  return 0;
}
