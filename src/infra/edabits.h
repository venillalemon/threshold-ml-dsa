// src/infra/edabits.h — demo edaBits, sampled and shared directly by the dealer.
//
// The FakeDealer samples each value ALREADY in the required range (no rejection
// sampling) and deals all needed representations of the SAME value at once:
//   * F_q share            (deal_fq)
//   * widened ring share   (deal_ring)
//   * authenticated GMW/WRK Boolean bits (deal_bits), consistent under the same
//     per-party Delta the backend uses.
// So an edaBit costs zero communication and zero GMW here; only the actual
// signing circuits (C_ad, C_prod, C_post) run real MPC. This is the paper's
// "KeyGen / preprocessing supplies the correlated setup state" assumption,
// realized as a trusted-dealer stand-in (the same cheat class as the dealer
// knowing every party's authentication key). A real edaBit protocol produces
// these pairs without a dealer and without opening the value.
#ifndef MLDSA_EDABITS_H
#define MLDSA_EDABITS_H

#include "backend.h"
#include "spdz.h"

#include <cstdint>
#include <vector>

namespace mldsa {

// Authenticated arithmetic and Boolean representations of the same value.
template <int nP, int BitWidth, bool IsUnsigned> struct SharePair {
  static_assert(BitWidth > 0, "SharePair: BitWidth must be positive");
  static constexpr int bit_width = BitWidth;
  static constexpr bool is_unsigned = IsUnsigned;

  std::vector<FqShare<nP>> fq_share;     // pairwise-authenticated F_q share
  std::vector<RingShare<nP>> ring_share; // same value over Z_{2^RING_W}
  emp::wrk::AuthShareVec<nP> two_share;  // GMW/WRK authenticated Boolean bits, LSB first
};

// Ring eDaBit (paper: eDaBit_{2^k,k}): a widened R^ in Z_{2^RING_W} and the low
// KB bits of R^ as authenticated Boolean shares.
template <int nP> struct RingEdabits {
  std::vector<RingShare<nP>> ring_share;
  emp::wrk::AuthShareVec<nP> two_share; // KB bits per value, little-endian
};

// y edaBit: y uniform in (-2^(W-1), 2^(W-1)] (ML-DSA's left-open interval).
// Boolean half is the two's-complement W-bit encoding of v = y - 1; arithmetic
// halves are y itself over F_q and the ring.
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, false> y_edabits(Backend<nP>& /*bk*/, int /*party*/,
                                                FakeDealer<nP>& dealer, int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "y_edabits: BitWidth in [1,31]");
  SharePair<nP, BitWidth, false> out;
  out.fq_share.resize((size_t)count);
  out.ring_share.resize((size_t)count);
  out.two_share.resize((size_t)BitWidth * count);
  for (int i = 0; i < count; ++i) {
    const uint64_t vu = dealer.draw_uniform(uint64_t{1} << BitWidth); // v = s (two's complement)
    const int64_t s = (vu >> (BitWidth - 1)) ? (int64_t)vu - (int64_t{1} << BitWidth) : (int64_t)vu;
    const int64_t y = s + 1;
    const uint32_t yq = (uint32_t)(((y % Q) + Q) % Q);
    out.fq_share[(size_t)i] = dealer.deal_fq(yq);
    out.ring_share[(size_t)i] = dealer.deal_ring(ring_from_i64(y));
    auto bits = dealer.deal_bits(vu, BitWidth);
    for (int k = 0; k < BitWidth; ++k)
      out.two_share[(size_t)BitWidth * i + k] = bits[(size_t)k];
  }
  return out;
}

// F_q edaBit: R uniform in [0, Q), with its canonical W-bit Boolean encoding.
// (ring_share is left empty; the A2B mask only needs F_q + Boolean.)
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, true> Fq_edabits(Backend<nP>& /*bk*/, int /*party*/,
                                                FakeDealer<nP>& dealer, int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "Fq_edabits: BitWidth in [1,31]");
  SharePair<nP, BitWidth, true> out;
  out.fq_share.resize((size_t)count);
  out.two_share.resize((size_t)BitWidth * count);
  for (int i = 0; i < count; ++i) {
    const uint32_t r = (uint32_t)dealer.draw_uniform((uint64_t)Q);
    out.fq_share[(size_t)i] = dealer.deal_fq(r);
    auto bits = dealer.deal_bits(r, BitWidth);
    for (int k = 0; k < BitWidth; ++k)
      out.two_share[(size_t)BitWidth * i + k] = bits[(size_t)k];
  }
  return out;
}

// Ring edaBit: widened R^ uniform in Z_{2^RING_W}; the low KB bits are the
// Boolean half.
template <int nP, int KB>
inline RingEdabits<nP> ring_edabits(Backend<nP>& /*bk*/, int /*party*/, FakeDealer<nP>& dealer,
                                    int count) {
  static_assert(KB > 0 && KB < 32, "ring_edabits: KB in [1,31]");
  constexpr uint64_t KMASK = (uint64_t{1} << KB) - 1;
  RingEdabits<nP> out;
  out.ring_share.resize((size_t)count);
  out.two_share.resize((size_t)KB * count);
  for (int i = 0; i < count; ++i) {
    const RingVal v = dealer.ring_rand();
    out.ring_share[(size_t)i] = dealer.deal_ring(v);
    auto bits = dealer.deal_bits(v.w[0] & KMASK, KB);
    for (int k = 0; k < KB; ++k)
      out.two_share[(size_t)KB * i + k] = bits[(size_t)k];
  }
  return out;
}

} // namespace mldsa
#endif // MLDSA_EDABITS_H
