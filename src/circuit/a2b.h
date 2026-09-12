// src/circuit/a2b.h — F_A2B: convert an authenticated arithmetic sharing <x>_q into
// Boolean circuit values <x>_2 inside the AG session, by mask-and-open:
//
//   1. draw edaBit masks (<r>_q, <r>_2) with r uniform in [0,Q)   (Fq_edabits)
//   2. open c = x + r mod q publicly, checked pairwise (BDOZ) in one flight
//   3. inside the circuit compute x = c - r mod q                 (sub_modq)
//
// c = x + r is uniform in F_q independently of x (r is a one-time pad), so the
// open leaks nothing about x; the MACCheck is what makes the public c
// malicious-secure — a tampered share flips c and the check aborts (see
// TAMPER_C). Output wires are unsigned 23-bit representatives in [0, Q),
// ready to feed decompose<PARAM>().
#ifndef MLDSA_A2B_H
#define MLDSA_A2B_H

#include "circuit.h" // fa, pack_bits/unpack_bits circuit primitives
#include "edabits.h" // Fq_edabits, SharePair
#include "spdz.h"      // FqShare, checked opens
#include <emp-ag/emp-ag.h>

#include <cstdint>
#include <random>
#include <vector>

namespace mldsa {

// ---------------------------------------------------------------------------
// A2B subtraction: w = (c - r) mod Q, both in [0, Q).  c is the SPDZ masked-open
// (PUBLIC), r is the adopted ⟨r⟩₂.  Reconstructs w INSIDE the circuit; feed w
// straight into decompose<PARAM>(ctx, w).  Reuses fa / pack_bits / unpack_bits
// and the "x - y == x + ~y + 1" ripple-subtract idiom (see decompose_g88).
//
// Contract: c, r, w all in [0, Q), Q < 2^L, little-endian L-bit. We work mod Q
// end-to-end (not over the integers), so no slack bits are needed.  Two L-bit
// ripple chains; c and Q are public constants, so only the carry chains cost AND.
// ---------------------------------------------------------------------------
template <class Ctx>
inline UInt_T<Ctx, L> sub_modq(Ctx& ctx, const UInt_T<Ctx, L>& c, // public
                               const UInt_T<Ctx, L>& r) {         // shared
  Bit_T<Ctx> cb[L], rb[L];
  unpack_bits(c, cb);
  unpack_bits(r, rb);

  // ---- Step A: subtract via add.  s = c + (~r) + 1  ==  c - r + 2^L. ----
  // The trailing +1 is the two's-complement "form -r" carry-in (ALWAYS on, not
  // a fixup). The carry OUT of bit 22 is the sign flag:
  //     cout == 1  <=>  c >= r   (result already in [0,Q), no fixup)
  //     cout == 0  <=>  c <  r   (result is negative, owe one +Q)
  Bit_T<Ctx> d[L], cout = Bit_T<Ctx>::constant(ctx, true); // carry-in = 1
  for (int i = 0; i < L; ++i)
    fa(d[i], cout, cb[i], !rb[i], cout); // d = (c - r) mod 2^L
  // (aliasing: s=d[i] distinct, cout==cin allowed — same rule fa documents.)

  // ---- Step B: conditional +Q when we borrowed (need = !cout). ----
  // Adding Q (NOT 1) is the mod-Q correction. Exactly one +Q suffices because
  // c-r in (-Q,0) => c-r+Q in (0,Q); the top carry is the 2^L that drops.
  const Bit_T<Ctx> need = !cout;
  Bit_T<Ctx> w[L], carry = Bit_T<Ctx>::constant(ctx, false); // carry-in = 0
  for (int i = 0; i < L; ++i) {
    // Q's bits are compile-time constants: where Q[i]=1 the addend is `need`
    // itself, where Q[i]=0 the full adder degenerates to a half adder. Folding
    // this saves one AND per bit over materialising need & Q[i].
    if (((Q >> i) & 1) != 0)
      fa(w[i], carry, d[i], need, carry);
    else
      ha(w[i], carry, d[i], carry);
  }
  // final carry dropped => w = (c - r) mod Q, in [0, Q).
  return pack_bits<L>(ctx, w);
}

template <int nP>
inline std::vector<emp::UInt_T<typename emp::AGMPCSession<nP>::ctx_t, L>>
a2b(emp::AGMPCSession<nP>& sess, int party, FakeDealer<nP>& dealer,
    const std::vector<FqShare<nP>>& x_share) {
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using U23 = emp::UInt_T<Ctx, L>;
  const int count = (int)x_share.size();

  // 1. one edaBit mask per coefficient.
  SharePair<nP, L, true> r = Fq_edabits<nP, L>(sess, party, dealer, count);

  // 2. open c = x + r with an SPDZ MACCheck.
  std::vector<FqShare<nP>> c_share((size_t)count);
  for (int i = 0; i < count; ++i)
    c_share[(size_t)i] = x_share[(size_t)i] + r.fq_share[(size_t)i];
#ifdef TAMPER_C
  if (party == 1 && count > 0)
    c_share[0].val ^= 1;
#endif
  const std::vector<uint32_t> c = open_fq_checked(sess.io(), party, dealer.my_alpha_f, c_share);

  // 3. x = c - r mod q inside the circuit.
  std::vector<U23> out;
  out.reserve((size_t)count);
  for (int i = 0; i < count; ++i) {
    U23 rr = sess.template adopt_authenticated_input<U23>(&r.two_share[(size_t)L * i]);
    U23 c_public = sess.template input<U23>(emp::PUBLIC, (uint64_t)c[i]);
    out.push_back(sub_modq(sess.ctx(), c_public, rr));
  }
  return out;
}

} // namespace mldsa
#endif // MLDSA_A2B_H
