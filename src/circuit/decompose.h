// src/circuit/decompose.h — Decompose for ML-DSA (FIPS 204) as a boolean circuit.
//
// The HighBits half and its gate count are ported verbatim from the emp-sh2pc
// version in files/mldsa_highbits.h; LowBits is new (it falls out of the same
// intermediates — see above decompose_g32). Only the API changed from the
// original: emp-sh2pc's `emp::Bit` rides a
// process-global backend, while emp-tool's `Bit_T<Ctx>` carries its context
// explicitly. So every function here takes a `Ctx&` and is templated on the
// BooleanContext — which is what lets the SAME code run on ClearCtx
// (plaintext), CountCtx (gate counting), and AGMPCSession<nP>::ctx_t (n-party
// malicious MPC).
//
//   decompose_g32()  gamma2 = (q-1)/32   65 / 87   w1 4 bits + w0 19   43 AND
//   decompose_g88()  gamma2 = (q-1)/88   44        w1 6 bits + w0 18  100 AND
//
// There is no separate HighBits: it is these with w0 dropped, which DCE takes
// back to 21 / 69 AND. See the comment above decompose_g32.
//
// ---------------------------------------------------------------------------
// Why it is this cheap: FACTOR the divisor instead of multiplying by its
// reciprocal.
//
//   alpha = 2*gamma2 = 523776 = 2^9 * 1023 = 512 * (2^10 - 1)
//
//   * the 2^9 factor  -> dividing by 512 is just "take bits 9 and up".  0 gates.
//   * the 1023 factor -> 1023 is a MERSENNE number, and 1024 = 1023 + 1, so
//     1024 == 1 (mod 1023) and the base-1024 digits can simply be summed.
//     This is casting out nines (1/9 in decimal), one fold, exact.
//
// Derivation for gamma2 = (q-1)/32, with r1 = floor((r + alpha/2 - 1)/alpha):
//
//   alpha/2 - 1 = 261887 = 511*512 + 255
//   split r = rho*512 + L,  L = r[0..8] < 512,  rho = r[9..22]
//
//     r + alpha/2 - 1 = (rho + 511)*512 + (L + 255)
//     r1 = floor( (rho + 511 + (L+255)/512) / 1023 )
//        = floor( (rho + 511 + g) / 1023 )          g = [L >= 257]
//
//   split rho = R_hi*1024 + R_lo,  R_lo = r[9..18], R_hi = r[19..22]
//   since 1024 = 1023 + 1:   rho = R_hi*1023 + (R_hi + R_lo)
//
//     r1 = R_hi + floor( (R_hi + R_lo + g + 511) / 1023 )
//        = R_hi + h                                 h = [n >= 512], n = R_lo+R_hi+g
//
//   n <= 1039 so n + 511 < 2*1023: the quotient is a SINGLE BIT. The answer is
//   literally bits 19..22 of r plus one correction bit.
//
// The rounding offset never becomes an addition: alpha/2 - 1 splits into 255
// and 511, which turn into the two comparison thresholds [L >= 257] and
// [n >= 512] — and a threshold is just a choice of which wire to look at.
//
// The FIPS 204 special case (r - r0 == q-1) is exactly r1 == 16. Since 16 is a
// power of two, "mod 16" drops it by never computing the carry out of bit 3 —
// the special case costs negative gates.
#ifndef MLDSA_DECOMPOSE_H
#define MLDSA_DECOMPOSE_H

#include "circuit.h" // fa/ha/or_tree/pack_bits circuit primitives
#include "ref.h"     // parameters + plaintext reference algorithms

#include <cstdint>
#include <emp-tool/emp-tool.h>
#include <vector>

namespace mldsa {

// Circuit output widths of Decompose (all other parameters live in ref.h).
constexpr int OW1 = (PARAM == MLDSA_44) ? 6 : 4;
constexpr int OW0 = (PARAM == MLDSA_44) ? 18 : 19;

// ---------------------------------------------------------------------------
// Decompose = HighBits *and* LowBits, for PrepSign: w1 goes out in the clear,
// w0 stays shared for the rejection test.
//
//   decompose_g32()  65 / 87   w1 4 bits, w0 19 bits two's complement   43 AND
//   decompose_g88()  44        w1 6 bits, w0 18 bits two's complement  100 AND
//
// There is no separate HighBits entry point: w1 alone is this circuit with w0
// dropped, and output-rooted DCE collapses it back to 21 / 69 AND. The AG
// backend garbles only live gates, so "call decompose, keep w1" costs exactly
// what a dedicated HighBits would.
//
// The obvious route for w0 is r0 = r - alpha*r1, a constant multiply plus a
// L-bit subtract. Don't: r0 falls out of the intermediates w1 ALREADY has.
//
// For gamma2 = (q-1)/32, in the notation of the derivation at the top of this
// file (rho = r[9..22], L = r[0..8], n = R_lo + R_hi + g, r1 = R_hi + h):
//
//   rho - 1023*r1 = R_hi*1023 + (R_hi + R_lo) - 1023*(R_hi + h)
//                 = (n - g) - 1023*h
//   r0 = r - alpha*r1 = 512*(rho - 1023*r1) + L
//                     = 512*(n - g - 1023*h) + L
//
// So r0 is m := n - g - 1023h in the high half and L verbatim in the low half —
// and since only 10 bits of m survive (r0 fits 19 bits two's complement),
// -1023h == +h (mod 1024) and the 1023 vanishes:  m = n - g + h  (mod 1024).
// Three +-1 corrections on a value already computed. No multiply, no subtract.
//
// The FIPS special case pays here, unlike in HighBits. When r1 would be 16,
// FIPS sets r1 := 0 AND r0 := r0 - 1 (keeping r == r1*alpha + r0 mod q). The
// flag is the carry HighBits threw away — 1 AND to keep it — and the decrement
// is a 9-bit borrow chain through L.
//
// Its cost is also why the flag must NOT be revealed to save those gates: it is
// public w1 == 0 that leaves w in two buckets, and the flag says which — i.e.
// it hands out sign(w0), which is exactly what the rejection test hides.
//
// ML-DSA-44 gets no such gift: alpha = 2^11 * 93 and the fold multiplies by 11,
// so backing r1 out costs a real 93*r1 (as 3*r1 shifted and subtracted, then
// truncated to the 7 bits that survive) — 29 AND on top of the 71.
// ---------------------------------------------------------------------------
template <class Ctx>
inline void decompose_g32(Ctx& ctx, Bit_T<Ctx>* w0, Bit_T<Ctx>* w1, const Bit_T<Ctx>* r) {
  (void)ctx;
  // stage 1: g = [L >= 257] = r[8] & (r[0] | ... | r[7])              8 AND
  //          L = r[0..8] < 512; L >= 257 iff bit 8 set and the low 8 nonzero.
  Bit_T<Ctx> g = r[8] & or_tree(r, 8);

  // stage 2: n = R_lo + R_hi + g = r[9..18] + r[19..22] + g, 11 bits  10 AND
  //          bits 0..3 see all three operands (full adders);
  //          bits 4..9 have only R_lo left plus the carry (half adders).
  Bit_T<Ctx> n[11], c;
  fa(n[0], c, r[9], r[19], g);
  for (int i = 1; i < 4; i++)
    fa(n[i], c, r[9 + i], r[19 + i], c);
  for (int i = 4; i < 10; i++)
    ha(n[i], c, r[9 + i], c);
  n[10] = c;

  // stage 3: h = [n >= 512] = n[9] | n[10].  n <= 1039, so bit 10 set means
  //          n in [1024,1039] and bit 9 is clear: the two are disjoint and
  //          the OR degenerates to a free XOR.                        0 AND
  Bit_T<Ctx> h = n[9] ^ n[10];

  // stage 4: w1 = (R_hi + h) mod 16, a 4-bit incrementer whose top carry is
  //          NOT dropped: that carry is "r1 == 16", the FIPS special case.
  //          Dropping it from w1 is the mod 16; keeping it is w0's job. 4 AND
  Bit_T<Ctx> k = h;
  for (int i = 0; i < 4; i++)
    ha(w1[i], k, r[19 + i], k);
  const Bit_T<Ctx> sp = k;

  // stage 5: low half, w0[0..8] = L - sp. A borrow chain, not an ha: the
  //          carry is `bo & !L[i]`, not `L[i] & bo`. bo is the borrow that
  //          then propagates into m.                                    9 AND
  Bit_T<Ctx> bo = sp;
  for (int i = 0; i < 9; i++) {
    w0[i] = r[i] ^ bo;
    bo = bo & !r[i];
  }

  // stage 6: s = h - g - bo in [-2,1], as a 2-bit two's complement (s1,s0)
  //          so it can be sign-extended into the 10-bit add below.      3 AND
  //          s1 = [s < 0] = h ? (g & bo) : (g | bo)
  Bit_T<Ctx> gb_and = g & bo, gb_or = g | bo;
  Bit_T<Ctx> s0 = h ^ g ^ bo;
  Bit_T<Ctx> s1 = gb_or ^ (h & (gb_or ^ gb_and));

  // stage 7: high half, w0[9..18] = m = (n + s) mod 2^10. The carry out of
  //          the top bit is dropped, which is what makes the mod free.  9 AND
  ha(w0[9], c, n[0], s0);
  for (int i = 1; i < 9; i++)
    fa(w0[9 + i], c, n[i], s1, c);
  w0[18] = n[9] ^ s1 ^ c;
}

// ---------------------------------------------------------------------------
// ML-DSA-44 :  gamma2 = (q-1)/88 = 95232,  alpha = 190464 = 2^11 * 93
//
// 93 is not Mersenne, so the fold does not apply directly. But 11*93 = 1023,
// so  floor(t/93) = floor(11t/1023)  and one pre-multiply by 11 restores it.
// alpha/2 - 1 = 95231 = 46*2048 + 1023.
//
//   in : r[0..22] little-endian    out: w1[0..5], w0[0..17]      100 AND
//                                       (w1 alone survives DCE at 71)
// ---------------------------------------------------------------------------
template <class Ctx>
inline void decompose_g88(Ctx& ctx, Bit_T<Ctx>* w0, Bit_T<Ctx>* w1, const Bit_T<Ctx>* r) {
  // stage 1: g = [L >= 1025] = r[10] & (r[0] | ... | r[9])           10 AND
  Bit_T<Ctx> g = r[10] & or_tree(r, 10);

  // stage 2: t = rho + 46 + g, rho = r[11..22], 13 bits              12 AND
  //          46 folds into the adder bit by bit — the carry is `x|c` where
  //          the constant bit is 1 and `x&c` where it is 0, so this is not
  //          an ha. g rides in as the carry-in for free.
  const int K = 46;
  Bit_T<Ctx> t[13], c = g;
  for (int i = 0; i < 12; i++) {
    const Bit_T<Ctx>& x = r[11 + i];
    if ((K >> i) & 1) {
      t[i] = !(x ^ c);
      c = x | c;
    } else {
      t[i] = (x ^ c);
      c = x & c;
    }
  }
  t[12] = c;

  // stage 3: u = 11t = 3t + (t<<3).  Two sequential ripple adds beat a
  //          CSA + CPA here (26 vs 28) because the operands are sparse. 26 AND
  Bit_T<Ctx> p[14]; // p = t + (t<<1) = 3t   13 AND
  p[0] = t[0];
  ha(p[1], c, t[1], t[0]);
  for (int i = 2; i < 13; i++)
    fa(p[i], c, t[i], t[i - 1], c);
  ha(p[13], c, t[12], c);

  Bit_T<Ctx> u[16]; // u = p + (t<<3)        13 AND
  for (int i = 0; i < 3; i++)
    u[i] = p[i];
  ha(u[3], c, p[3], t[0]);
  for (int i = 4; i < 14; i++)
    fa(u[i], c, p[i], t[i - 3], c);
  ha(u[14], c, t[11], c);
  ha(u[15], c, t[12], c);

  const Bit_T<Ctx>* b = u;      // u[0..9]   the low digit
  const Bit_T<Ctx>* a = u + 10; // u[10..15] the high digit, <= 44

  // stage 4: h = [a + b >= 1023] = carry out of bit 10 of (a + b + 1).
  //          Only the carry chain is needed; every sum bit is discarded. 10 AND
  Bit_T<Ctx> dump;
  fa(dump, c, b[0], a[0], Bit_T<Ctx>::constant(ctx, true));
  for (int i = 1; i < 6; i++)
    fa(dump, c, b[i], a[i], c);
  for (int i = 6; i < 10; i++)
    ha(dump, c, b[i], c);
  Bit_T<Ctx> h = c;

  // stage 5: r1 = a + h, 6 bits; the top carry is unused. This is the RAW
  //          quotient — LowBits needs it un-masked, w1 needs it masked. 5 AND
  Bit_T<Ctx> r1[6], k = h;
  for (int i = 0; i < 5; i++)
    ha(r1[i], k, a[i], k);
  r1[5] = a[5] ^ k;

  // stage 6: r1 == 44 (0b101100) -> 0.  Unlike the 65/87 case this costs real
  //          gates, because 44 is not a power of two.                   8 AND
  Bit_T<Ctx> f = ((r1[5] & !r1[4]) & (r1[3] & r1[2])) & (!r1[1] & !r1[0]);
  Bit_T<Ctx> nf = !f;
  for (int i = 0; i < 6; i++)
    w1[i] = r1[i];
  w1[2] = r1[2] & nf; // bits 0,1,4 are already 0 when r1 == 44,
  w1[3] = r1[3] & nf; // so only the set bits of 44 need masking
  w1[5] = r1[5] & nf;

  // ---- LowBits: w0 = r - alpha*r1 - f,  alpha = 2^11 * 93 ---------- 29 AND
  // Only 18 bits of w0 survive (|w0| <= gamma2 = 95232), so above bit 11 only
  // 7 bits of 93*r1 are ever needed — everything below is truncated away.

  Bit_T<Ctx> v[8]; // v = 3*r1                6 AND
  v[0] = r1[0];
  ha(v[1], c, r1[1], r1[0]);
  for (int i = 2; i < 6; i++)
    fa(v[i], c, r1[i], r1[i - 1], c);
  ha(v[6], c, r1[5], c);
  v[7] = c;

  Bit_T<Ctx> N[7];      // N = -v mod 2^7          5 AND
  N[0] = v[0];          // copy up to the lowest set
  Bit_T<Ctx> kk = v[0]; // bit, invert above it
  for (int i = 1; i < 7; i++) {
    N[i] = v[i] ^ kk;
    if (i < 6)
      kk = v[i] | kk;
  }

  Bit_T<Ctx> D[7]; // D = 93*r1 = 31*v = 32v-v  1 AND
  for (int i = 0; i < 5; i++)
    D[i] = N[i];           // mod 2^7, so 32v is just
  ha(D[5], c, N[5], v[0]); // v[0],v[1] at bits 5,6
  D[6] = N[6] ^ v[1] ^ c;

  Bit_T<Ctx> bo = f; // w0[0..10] = L - f      11 AND
  for (int i = 0; i < 11; i++) {
    w0[i] = r[i] ^ bo;
    bo = bo & !r[i];
  }

  // w0[11..17] = (r>>11) - D - bo.  As -D-bo == ~D + !bo, the borrow rides in
  // as the carry-in and it is one plain ripple add; the carry out is the
  // truncation and is dropped.                                          6 AND
  fa(w0[11], c, r[11], !D[0], !bo);
  for (int i = 1; i < 6; i++)
    fa(w0[11 + i], c, r[11 + i], !D[i], c);
  w0[17] = r[17] ^ !D[6] ^ c;
}

// ---------------------------------------------------------------------------
// Typed wrappers — what a session actually calls.
// UInt_T<Ctx,L> in, {w1, w0} out. The pack/unpack is pure rewiring, so these
// cost exactly the same 43 / 100 AND.
//
// For HighBits alone, call these and drop w0: DCE takes it back to 21 / 69.
// Under a chunked session "drop" means let the value leave C++ scope before the
// flush — reveal keeps every wire still in scope (live_pending_ids), so a w0
// you are still holding is a w0 the backend still garbles.
// ---------------------------------------------------------------------------

// w0 is a TWO'S COMPLEMENT value in a UInt_T — the natural encoding for the
// rejection test, since |w0| < gamma2 <= 2^(OW0-1) always holds and additions
// downstream are the same gates either way.
template <class Ctx, int OW1, int OW0> struct Decomposed {
  UInt_T<Ctx, OW1> w1; // revealed
  UInt_T<Ctx, OW0> w0; // stays shared
};

template <class Ctx>
inline Decomposed<Ctx, 4, 19> decompose_g32(Ctx& ctx, const UInt_T<Ctx, L>& r) {
  Bit_T<Ctx> rb[L], w1b[4], w0b[19];
  unpack_bits(r, rb);
  decompose_g32(ctx, w0b, w1b, rb);
  return {pack_bits<4>(ctx, w1b), pack_bits<19>(ctx, w0b)};
}

template <class Ctx>
inline Decomposed<Ctx, 6, 18> decompose_g88(Ctx& ctx, const UInt_T<Ctx, L>& r) {
  Bit_T<Ctx> rb[L], w1b[6], w0b[18];
  unpack_bits(r, rb);
  decompose_g88(ctx, w0b, w1b, rb);
  return {pack_bits<6>(ctx, w1b), pack_bits<18>(ctx, w0b)};
}

// Parameter-set dispatch. Must be a template so `if constexpr` actually
// discards the unused branch — the two paths have different return widths.
template <ParamSet P, class Ctx> inline auto decompose(Ctx& ctx, const UInt_T<Ctx, L>& r) {
  if constexpr (P == MLDSA_44)
    return decompose_g88(ctx, r);
  else
    return decompose_g32(ctx, r);
}

} // namespace mldsa
#endif // MLDSA_DECOMPOSE_H
