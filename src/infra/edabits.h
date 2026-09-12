// src/infra/edabits.h — demo edaBits for the PrepSign prototype.
//
// Returns authenticated arithmetic and Boolean representations. The unsigned
// and signed edaBit APIs represent the same value in both domains; y_edabits
// uses the documented affine encoding two_share=y-1, q_share=y.
//
// This is deliberately the same CHEAT that used to live in main.cpp: draw the
// Boolean shares, XOR-open r, then dealer-split that public r in F_Q. A real
// implementation must produce both representations without opening r.
#ifndef MLDSA_EDABITS_H
#define MLDSA_EDABITS_H

#include "spdz.h"
#include <emp-ag/emp-ag.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

namespace mldsa {

// Authenticated arithmetic and Boolean representations. They match exactly
// for unsigned_edabits, signed_edabits, Fq_edabits and rand. y_edabits is the
// one documented affine-encoding exception. BitWidth and IsUnsigned are part
// of the type, so consumers cannot silently mix incompatible encodings.
template <int nP, int BitWidth, bool IsUnsigned> struct SharePair {
  static_assert(BitWidth > 0, "SharePair: BitWidth must be positive");
  static constexpr int bit_width = BitWidth;
  static constexpr bool is_unsigned = IsUnsigned;

  std::vector<FqShare<nP>> fq_share;     // pairwise-authenticated F_q share
  std::vector<RingShare<nP>> ring_share; // same value in Z_{2^RING_W} (two's complement if signed)
  emp::ag::AShareBundleVec<nP> two_share;
};

// Ring eDaBit (paper: eDaBit_{2^k,k}): a uniform widened R^ in Z_{2^RING_W}
// as a ring share, and the low KB bits of R^ as Boolean shares.
template <int nP> struct RingEdabits {
  std::vector<RingShare<nP>> ring_share;
  emp::ag::AShareBundleVec<nP> two_share; // KB bits per value, little-endian
};

template <int nP, int BitWidth, bool IsUnsigned>
inline SharePair<nP, BitWidth, IsUnsigned> edabits(emp::AGMPCSession<nP>& sess, int party,
                                                   FakeDealer<nP>& dealer, int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "edabits: BitWidth must be in [1,31]");
  emp::expecting(count >= 0, "edabits: count must be non-negative");

  SharePair<nP, BitWidth, IsUnsigned> out;
  if (count == 0)
    return out;
  sess.protocol().fpre->abit.draw(BitWidth * count, out.two_share);

  // Pack this party's local XOR-bit shares for each coefficient. The local
  // share bit is encoded in bit 0 of every peer MAC, so slot 0 is sufficient.
  std::vector<uint32_t> xor_r_share((size_t)count);
  for (int i = 0; i < count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < BitWidth; ++k)
      v |= (uint32_t)emp::getLSB(out.two_share[BitWidth * i + k].mac(0)) << k;
    xor_r_share[i] = v;
  }

  // DEMO CHEAT: reveal r only to manufacture a matching arithmetic share.
  // The returned object does not expose this public intermediate.
  const std::vector<uint32_t> r = open_xor(sess.io(), party, xor_r_share);
  out.fq_share.resize((size_t)count);
  out.ring_share.resize((size_t)count);
  constexpr uint64_t VALUE_RANGE = uint64_t{1} << BitWidth;
  constexpr uint32_t RANGE_MOD_Q = (uint32_t)(VALUE_RANGE % Q);
  constexpr uint32_t SIGN_BIT = uint32_t{1} << (BitWidth - 1);
  for (int i = 0; i < count; ++i) {
    uint32_t value = r[i] % Q;
    RingVal ring_value = ring_from_u64(r[i]);
    if constexpr (!IsUnsigned)
      if (r[i] & SIGN_BIT) {
        value = fq_sub(value, RANGE_MOD_Q);
        ring_value = ring_sub(ring_value, ring_from_u64(VALUE_RANGE)); // sign-extend
      }
    out.fq_share[i] = dealer.deal_fq(value);
    out.ring_share[i] = dealer.deal_ring(ring_value);
  }

  return out;
}

// Ring eDaBits: draw RING_K Boolean bits per value, XOR-open them (same DEMO
// cheat as above), and let the dealer split a widened R^ whose low RING_K
// bits equal the opened R and whose high bits come from the shared stream.
template <int nP, int KB>
inline RingEdabits<nP> ring_edabits(emp::AGMPCSession<nP>& sess, int party,
                                    FakeDealer<nP>& dealer, int count) {
  static_assert(KB > 0 && KB < 32, "ring_edabits: KB must be in [1,31]");
  emp::expecting(count >= 0, "ring_edabits: count must be non-negative");
  RingEdabits<nP> out;
  if (count == 0)
    return out;
  constexpr uint64_t KMASK = (uint64_t{1} << KB) - 1;
  sess.protocol().fpre->abit.draw(KB * count, out.two_share);
  std::vector<uint32_t> xor_r_share((size_t)count);
  for (int i = 0; i < count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < KB; ++k)
      v |= (uint32_t)emp::getLSB(out.two_share[(size_t)KB * i + k].mac(0)) << k;
    xor_r_share[(size_t)i] = v;
  }
  const std::vector<uint32_t> r = open_xor(sess.io(), party, xor_r_share);
  out.ring_share.resize((size_t)count);
  for (int i = 0; i < count; ++i) {
    RingVal v = dealer.ring_rand(); // uniform widened lift R^
    v.w[0] = (v.w[0] & ~KMASK) | (uint64_t)r[(size_t)i]; // low KB bits are the Boolean half
    out.ring_share[(size_t)i] = dealer.deal_ring(v);
  }
  return out;
}

// Unsigned edaBits over [0, 2^BitWidth).
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, true> unsigned_edabits(emp::AGMPCSession<nP>& sess, int party,
                                                      FakeDealer<nP>& dealer, int count) {
  return edabits<nP, BitWidth, true>(sess, party, dealer, count);
}

// Signed edaBits over [-2^(BitWidth-1), 2^(BitWidth-1)-1].
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, false> signed_edabits(emp::AGMPCSession<nP>& sess, int party,
                                                     FakeDealer<nP>& dealer, int count) {
  return edabits<nP, BitWidth, false>(sess, party, dealer, count);
}

// If signed_edabits returns s in
// [-2^(BitWidth-1), 2^(BitWidth-1)-1], then y=s+1 lies in
// (-2^(BitWidth-1), 2^(BitWidth-1)]. The Boolean output deliberately remains
// the two's-complement sharing of s=y-1; the arithmetic output is a sharing of
// y itself. A later Boolean circuit can therefore compute with y as s+1 while
// retaining the exact left-open/right-closed ML-DSA sampling interval.
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, false> y_edabits(emp::AGMPCSession<nP>& sess, int party,
                                                FakeDealer<nP>& dealer, int count) {
  SharePair<nP, BitWidth, false> out = signed_edabits<nP, BitWidth>(sess, party, dealer, count);
  // y = s + 1: add a public constant the SPDZ way — one party adds it to the
  // value share, EVERY party adds alpha_p*1 to its MAC share (no full alpha).
  for (FqShare<nP>& fs : out.fq_share)
    fs = add_public_fq(fs, 1, party, dealer.my_alpha_f);
  for (RingShare<nP>& rs : out.ring_share)
    rs = add_public_ring(rs, 1, party, dealer.my_alpha_r);
  return out;
}

// edaBits of a uniform element of F_q: draw unsigned BitWidth-bit edaBits and
// reject those >= Q (oversampled so the batch suffices except with probability
// 2^-lambda).
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, true> Fq_edabits(emp::AGMPCSession<nP>& sess, int party,
                                                FakeDealer<nP>& dealer, int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "Fq_edabits: BitWidth must be in [1,31]");
  emp::expecting(count >= 0, "Fq_edabits: count must be non-negative");
  if (count == 0)
    return {};

  constexpr uint64_t VALUE_RANGE = uint64_t{1} << BitWidth;
  if constexpr (VALUE_RANGE <= (uint64_t)Q)
    return unsigned_edabits<nP, BitWidth>(sess, party, dealer, count);

  const int sample_count = (int)std::ceil(
      (double)count * VALUE_RANGE * (1.0 + std::sqrt(2 * lambda * std::log(2.0) / count)) / Q);
  SharePair<nP, BitWidth, true> out =
      unsigned_edabits<nP, BitWidth>(sess, party, dealer, sample_count);

  std::vector<uint32_t> xor_r_share((size_t)sample_count);
  for (int i = 0; i < sample_count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < BitWidth; ++k)
      v |= (uint32_t)emp::getLSB(out.two_share[BitWidth * i + k].mac(0)) << k;
    xor_r_share[i] = v;
  }
  const std::vector<uint32_t> r = open_xor(sess.io(), party, xor_r_share);
  std::vector<int> valid_index((size_t)sample_count);
  std::iota(valid_index.begin(), valid_index.end(), 0);
  valid_index.erase(std::remove_if(valid_index.begin(), valid_index.end(),
                                   [&](int i) { return r[i] >= (uint32_t)Q; }),
                    valid_index.end());
  emp::expecting((int)valid_index.size() >= count,
                 "Fq_edabits: oversized batch contained too few values in Z_Q");

  for (int i = 0; i < count; ++i) {
    const int src = valid_index[i];
    out.fq_share[i] = out.fq_share[src];
    for (int k = 0; k < BitWidth; ++k)
      out.two_share[BitWidth * i + k] = out.two_share[BitWidth * src + k];
  }
  out.fq_share.resize((size_t)count);
  out.two_share.resize((size_t)BitWidth * count);

  return out;
}

} // namespace mldsa
#endif // MLDSA_EDABITS_H
