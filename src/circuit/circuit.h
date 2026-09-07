// src/circuit/circuit.h — Boolean circuit micro-components shared by decompose.h,
// a2b.h and sign.h.
//
// XOR / NOT are free under free-XOR; only & and | cost an AND gate. Everything
// is templated on the BooleanContext, so the SAME code runs on ClearCtx
// (plaintext), CountCtx (gate counting), and AGMPCSession<nP>::ctx_t (n-party
// malicious MPC).
//
// Ripple carry chains everywhere, on purpose. In garbled circuits the whole
// netlist ships in one shot, so AND *depth* is free and only AND *count* pays;
// a parallel-prefix adder trades 10 AND for 22 AND to buy depth nobody spends.
// (Retargeting to GMW would flip this — swap the chains for Kogge-Stone.)
#ifndef MLDSA_CIRCUIT_H
#define MLDSA_CIRCUIT_H

#include <cstdint>
#include <emp-tool/emp-tool.h>
#include <vector>

namespace mldsa {

using emp::Bit_T;
using emp::UInt_T;

// ---- adders ---------------------------------------------------------------
// ALIASING CONTRACT for fa/ha: `cout` MAY be the same object as `cin`
// — every ripple chain below is written as fa(sum[i], c, x, y, c), which is
// safe because the right-hand side is fully evaluated before the assignment,
// and t1/t2 have already snapshotted cin. But `s` must NOT alias cin (or a/b):
// writing s first would clobber cin before cout reads it.
//
// Reassigning a carry in a loop is not a copy that loses anything. Bit_T is a
// handle (a Ctx* and a wire id); the gate is recorded in the context when the
// operator runs, and rebinding the variable only changes what the NAME points
// at. Under ChunkRecorderCtx it also unpins the old id, which is exactly right:
// the refcount is the DCE root set, not reachability.

// 1-AND full adder.  s = a^b^cin,  cout = MAJ(a,b,cin).
template <class Ctx>
inline void fa(Bit_T<Ctx>& s, Bit_T<Ctx>& cout, const Bit_T<Ctx>& a, const Bit_T<Ctx>& b,
               const Bit_T<Ctx>& cin) {
  Bit_T<Ctx> t1 = a ^ cin, t2 = b ^ cin;
  s = t1 ^ b;
  cout = cin ^ (t1 & t2); // 1 AND
}

// 1-AND half adder (second operand implicitly 0).
template <class Ctx>
inline void ha(Bit_T<Ctx>& s, Bit_T<Ctx>& cout, const Bit_T<Ctx>& a, const Bit_T<Ctx>& cin) {
  s = a ^ cin;
  cout = a & cin; // 1 AND
}

// out = a + b + cin mod 2^W. out may alias a (accumulator use), not b.
template <int W, class Ctx>
inline void add_ripple(Ctx& ctx, Bit_T<Ctx>* out, const Bit_T<Ctx>* a, const Bit_T<Ctx>* b,
                       bool cin) {
  Bit_T<Ctx> c = Bit_T<Ctx>::constant(ctx, cin);
  for (int i = 0; i < W - 1; ++i)
    fa(out[i], c, a[i], b[i], c);
  out[W - 1] = a[W - 1] ^ b[W - 1] ^ c;
}

// out = a + b + cin mod 2^W with a circuit-bit carry-in. out may alias a, not b.
template <int W, class Ctx>
inline void add_ripple_cin(Ctx&, Bit_T<Ctx>* out, const Bit_T<Ctx>* a, const Bit_T<Ctx>* b,
                           const Bit_T<Ctx>& cin) {
  Bit_T<Ctx> c = cin;
  for (int i = 0; i < W - 1; ++i)
    fa(out[i], c, a[i], b[i], c);
  out[W - 1] = a[W - 1] ^ b[W - 1] ^ c;
}

// out = a + cin mod 2^W, cin a circuit bit (W-1 AND). out may alias a.
template <int W, class Ctx>
inline void add_bit_ripple(Ctx&, Bit_T<Ctx>* out, const Bit_T<Ctx>* a, const Bit_T<Ctx>& cin) {
  Bit_T<Ctx> c = cin;
  for (int i = 0; i < W - 1; ++i) {
    const Bit_T<Ctx> ai = a[i]; // snapshot: out may alias a
    out[i] = ai ^ c;
    c = ai & c;
  }
  out[W - 1] = a[W - 1] ^ c;
}

// out = a + k + cin mod 2^W, k a public constant. out may alias a.
template <int W, class Ctx>
inline void add_const_ripple(Ctx& ctx, Bit_T<Ctx>* out, const Bit_T<Ctx>* a, uint64_t k, bool cin) {
  Bit_T<Ctx> c = Bit_T<Ctx>::constant(ctx, cin);
  for (int i = 0; i < W; ++i) {
    const bool kb = (k >> i) & 1;
    const Bit_T<Ctx> ai = a[i]; // snapshot: out may alias a
    out[i] = kb ? !(ai ^ c) : (ai ^ c);
    if (i < W - 1)
      c = kb ? (ai | c) : (ai & c);
  }
}

// ---- comparisons ----------------------------------------------------------

// [v >= C] for signed W-bit v: sign bit of (v - C) at width W+1.
template <int W, class Ctx> inline Bit_T<Ctx> ge_const(Ctx& ctx, const Bit_T<Ctx>* v, int64_t C) {
  Bit_T<Ctx> ext[W + 1], t[W + 1];
  for (int i = 0; i < W; ++i)
    ext[i] = v[i];
  ext[W] = v[W - 1];
  add_const_ripple<W + 1>(ctx, t, ext, (uint64_t)(-C), false);
  return !t[W];
}

// [v >= C] for UNSIGNED W-bit v (C in [0, 2^W)): sign bit of (v - C) at width W+1.
template <int W, class Ctx>
inline Bit_T<Ctx> ge_const_u(Ctx& ctx, const Bit_T<Ctx>* v, uint64_t C) {
  Bit_T<Ctx> ext[W + 1], t[W + 1];
  for (int i = 0; i < W; ++i)
    ext[i] = v[i];
  ext[W] = Bit_T<Ctx>::constant(ctx, false);
  add_const_ripple<W + 1>(ctx, t, ext, (uint64_t)(-(int64_t)C), false);
  return !t[W];
}

// ---- reductions & rewiring ------------------------------------------------

// Balanced OR reduction: n-1 AND, depth ceil(log2 n).
template <class Ctx> inline Bit_T<Ctx> or_tree(const Bit_T<Ctx>* v, int n) {
  std::vector<Bit_T<Ctx>> cur(v, v + n);
  while (cur.size() > 1) {
    std::vector<Bit_T<Ctx>> nxt;
    size_t i = 0;
    for (; i + 1 < cur.size(); i += 2)
      nxt.push_back(cur[i] | cur[i + 1]);
    if (i < cur.size())
      nxt.push_back(cur[i]);
    cur.swap(nxt);
  }
  return cur[0];
}

// [v == C] for signed W-bit v, C in two's complement.
template <int W, class Ctx> inline Bit_T<Ctx> eq_const(Ctx& ctx, const Bit_T<Ctx>* v, int64_t C) {
  Bit_T<Ctx> diff[W];
  for (int i = 0; i < W; ++i)
    diff[i] = ((C >> i) & 1) ? !v[i] : v[i]; // v XOR C-bit: 1 where bits differ
  return !or_tree(diff, W);
}

// Bit array <-> typed value. Pure rewiring, 0 gates.
template <int N, class Ctx> inline UInt_T<Ctx, N> pack_bits(Ctx& ctx, const Bit_T<Ctx>* b) {
  typename Ctx::Wire w[N];
  for (int i = 0; i < N; ++i)
    b[i].pack_wires(w + i);
  return UInt_T<Ctx, N>::from_wires(ctx, w);
}
template <int N, class Ctx> inline void unpack_bits(const UInt_T<Ctx, N>& v, Bit_T<Ctx>* out) {
  for (int i = 0; i < N; ++i)
    out[i] = v[i];
}

} // namespace mldsa
#endif // MLDSA_CIRCUIT_H
