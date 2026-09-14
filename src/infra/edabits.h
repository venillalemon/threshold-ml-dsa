// src/infra/edabits.h — demo edaBits for the mixed-security PrepSign.
//
// Returns authenticated arithmetic and Boolean representations. The Boolean
// half is now an authenticated GMW/WRK share (emp::wrk::AuthShare) under the
// backend's single Delta, so it can be fed directly into an offline GMW circuit
// (C_pre) or installed as a fixed input of the WRK circuit (C_post). The
// arithmetic half is the FakeDealer BDOZ share, unchanged.
//
// This is deliberately the same CHEAT as before: draw the Boolean shares from
// GMW, XOR-open r, then dealer-split that public r in F_Q / the ring. A real
// implementation must produce both representations without opening r (e.g. a
// real edaBit protocol). Only the Boolean-share PROVENANCE changed (KRRW aBit
// pool -> authenticated GMW), not the cheat.
#ifndef MLDSA_EDABITS_H
#define MLDSA_EDABITS_H

#include "backend.h"
#include "spdz.h"

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
  emp::wrk::AuthShareVec<nP> two_share;  // GMW/WRK authenticated Boolean bits, LSB first
};

// Ring eDaBit (paper: eDaBit_{2^k,k}): a uniform widened R^ in Z_{2^RING_W}
// as a ring share, and the low KB bits of R^ as authenticated GMW Boolean bits.
template <int nP> struct RingEdabits {
  std::vector<RingShare<nP>> ring_share;
  emp::wrk::AuthShareVec<nP> two_share; // KB bits per value, little-endian
};

// The local XOR-share bit of an authenticated GMW/WRK Boolean share.
template <int nP> inline uint8_t local_bit(const emp::wrk::AuthShare<nP>& s) { return s.bit & 1; }

template <int nP, int BitWidth, bool IsUnsigned>
inline SharePair<nP, BitWidth, IsUnsigned> edabits(Backend<nP>& bk, int party,
                                                   FakeDealer<nP>& dealer, int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "edabits: BitWidth must be in [1,31]");
  emp::expecting(count >= 0, "edabits: count must be non-negative");

  SharePair<nP, BitWidth, IsUnsigned> out;
  if (count == 0)
    return out;
  out.two_share = bk.gmw().random((size_t)BitWidth * count);

  // Pack this party's local XOR-bit shares for each coefficient.
  std::vector<uint32_t> xor_r_share((size_t)count);
  for (int i = 0; i < count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < BitWidth; ++k)
      v |= (uint32_t)local_bit<nP>(out.two_share[BitWidth * i + k]) << k;
    xor_r_share[i] = v;
  }

  // DEMO CHEAT: reveal r only to manufacture a matching arithmetic share.
  const std::vector<uint32_t> r = open_xor(bk.io(), party, xor_r_share);
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

// Ring eDaBits: draw RING_K Boolean bits per value from GMW, XOR-open them (same
// DEMO cheat), and let the dealer split a widened R^ whose low RING_K bits equal
// the opened R and whose high bits come from the shared stream.
template <int nP, int KB>
inline RingEdabits<nP> ring_edabits(Backend<nP>& bk, int party, FakeDealer<nP>& dealer, int count) {
  static_assert(KB > 0 && KB < 32, "ring_edabits: KB must be in [1,31]");
  emp::expecting(count >= 0, "ring_edabits: count must be non-negative");
  RingEdabits<nP> out;
  if (count == 0)
    return out;
  constexpr uint64_t KMASK = (uint64_t{1} << KB) - 1;
  out.two_share = bk.gmw().random((size_t)KB * count);
  std::vector<uint32_t> xor_r_share((size_t)count);
  for (int i = 0; i < count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < KB; ++k)
      v |= (uint32_t)local_bit<nP>(out.two_share[(size_t)KB * i + k]) << k;
    xor_r_share[(size_t)i] = v;
  }
  const std::vector<uint32_t> r = open_xor(bk.io(), party, xor_r_share);
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
inline SharePair<nP, BitWidth, true> unsigned_edabits(Backend<nP>& bk, int party,
                                                      FakeDealer<nP>& dealer, int count) {
  return edabits<nP, BitWidth, true>(bk, party, dealer, count);
}

// Signed edaBits over [-2^(BitWidth-1), 2^(BitWidth-1)-1].
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, false> signed_edabits(Backend<nP>& bk, int party,
                                                     FakeDealer<nP>& dealer, int count) {
  return edabits<nP, BitWidth, false>(bk, party, dealer, count);
}

// If signed_edabits returns s in [-2^(BitWidth-1), 2^(BitWidth-1)-1], then
// y=s+1 lies in (-2^(BitWidth-1), 2^(BitWidth-1)]. The Boolean output stays the
// two's-complement sharing of s=y-1; the arithmetic output is a sharing of y.
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, false> y_edabits(Backend<nP>& bk, int party,
                                                FakeDealer<nP>& dealer, int count) {
  SharePair<nP, BitWidth, false> out = signed_edabits<nP, BitWidth>(bk, party, dealer, count);
  // y = s + 1: add a public constant the SPDZ way — one party adds it to the
  // value share, EVERY party adds alpha_p*1 to its MAC share (no full alpha).
  for (FqShare<nP>& fs : out.fq_share)
    fs = add_public_fq(fs, 1, party, dealer.my_alpha_f);
  for (RingShare<nP>& rs : out.ring_share)
    rs = add_public_ring(rs, 1, party, dealer.my_alpha_r);
  return out;
}

// edaBits of a uniform element of F_q: draw unsigned BitWidth-bit edaBits and
// reject those >= Q (oversampled so the batch suffices except w.p. 2^-lambda).
template <int nP, int BitWidth>
inline SharePair<nP, BitWidth, true> Fq_edabits(Backend<nP>& bk, int party, FakeDealer<nP>& dealer,
                                                int count) {
  static_assert(BitWidth > 0 && BitWidth < 32, "Fq_edabits: BitWidth must be in [1,31]");
  emp::expecting(count >= 0, "Fq_edabits: count must be non-negative");
  if (count == 0)
    return {};

  constexpr uint64_t VALUE_RANGE = uint64_t{1} << BitWidth;
  if constexpr (VALUE_RANGE <= (uint64_t)Q)
    return unsigned_edabits<nP, BitWidth>(bk, party, dealer, count);

  const int sample_count = (int)std::ceil(
      (double)count * VALUE_RANGE * (1.0 + std::sqrt(2 * lambda * std::log(2.0) / count)) / Q);
  SharePair<nP, BitWidth, true> out = unsigned_edabits<nP, BitWidth>(bk, party, dealer, sample_count);

  std::vector<uint32_t> xor_r_share((size_t)sample_count);
  for (int i = 0; i < sample_count; ++i) {
    uint32_t v = 0;
    for (int k = 0; k < BitWidth; ++k)
      v |= (uint32_t)local_bit<nP>(out.two_share[BitWidth * i + k]) << k;
    xor_r_share[i] = v;
  }
  const std::vector<uint32_t> r = open_xor(bk.io(), party, xor_r_share);
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
