// src/protocol/phase1_util.h — helpers to assemble offline (C_pre) circuits as
// emp-tool BooleanPrograms and run them through authenticated GMW.
//
// The per-coefficient circuit kernels (sub_modq, decompose, window_producer)
// are all templated on the BooleanContext, so they record straight into a
// RecordCtx. These helpers wrap consecutive input-wire ids as typed values and
// splice GMW public bits, so a whole Phase-1 stage compiles once and evaluates
// as one gmw.evaluate() call.
#ifndef MLDSA_PHASE1_UTIL_H
#define MLDSA_PHASE1_UTIL_H

#include "backend.h"
#include "circuit.h"

#include <emp-tool/circuits/frontend/frontend.h>
#include <emp-tool/ir/context/record.h>
#include <emp-tool/ir/program.h>

#include <cstdint>
#include <vector>

namespace mldsa {

// A slice of W consecutive input wires [base, base+W) as a typed value.
template <int W> inline emp::UInt_T<emp::RecordCtx, W> in_uint(emp::RecordCtx& ctx, uint32_t base) {
  emp::RecordCtx::Wire w[W];
  for (int i = 0; i < W; ++i)
    w[i] = base + (uint32_t)i;
  return emp::UInt_T<emp::RecordCtx, W>::from_wires(ctx, w);
}
template <class Ctx> inline emp::Bit_T<Ctx> wire_bit(Ctx& ctx, uint32_t id) {
  return emp::Bit_T<Ctx>(ctx, id);
}

// Append a value's W bits to a GMW-share input vector as public_bit shares.
template <int nP>
inline void push_public_word(Backend<nP>& bk, emp::wrk::ShareVec<nP>& in, uint64_t v, int W) {
  for (int k = 0; k < W; ++k)
    in.push_back(bk.gmw().public_bit((v >> k) & 1));
}

} // namespace mldsa

#endif // MLDSA_PHASE1_UTIL_H
