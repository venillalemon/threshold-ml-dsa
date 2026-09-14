// src/protocol/rand.h — F_smallnormpoly using daBits and authenticated GMW.
//
// For each output coefficient this returns matching authenticated shares:
//   q_share[i]                   = <r_i>_q over F_Q
//   ring_share[i]                = <r_i> over Z_{2^RING_W}
//   two_share[WIDTH*i + k]       = <bit_k(r_i)>_2 (GMW/WRK share), little-endian
// where WIDTH is 3 for eta=2 and 4 for eta=4, two's-complement encoded so
//   eta=2: b0 + 2*b1 - 4*b2
//   eta=4: b0 + 2*b1 + 4*b2 - 8*b3.
//
// The small-norm membership predicate is evaluated by GMW (the offline circuit
// domain): a couple of AND layers per candidate, then the reject bits open.
#ifndef MLDSA_RAND_H
#define MLDSA_RAND_H

#include "edabits.h"

#include <cstdint>
#include <cstdio>
#include <vector>

namespace mldsa {

// AND gates spent inside eDaBit/daBit rejection sampling. These are correlation
// generation, not protocol circuit work, so cost tables report them separately.
inline uint64_t g_edabit_ands = 0;

constexpr int64_t isqrt_floor(int64_t v) { // floor(sqrt(v)), constexpr-safe
  int64_t r = 0;
  while ((r + 1) * (r + 1) <= v)
    ++r;
  return r;
}

template <int nP, int ETA, int COUNT>
inline SharePair<nP, (ETA == 2 ? 3 : 4), false> rand_edabits(Backend<nP>& bk, int party,
                                                             FakeDealer<nP>& dealer) {
  static_assert(ETA == 2 || ETA == 4, "rand: ETA must be 2 or 4");
  static_assert(COUNT >= 0, "rand: COUNT must be non-negative");

  constexpr int WIDTH = ETA == 2 ? 3 : 4;
  // Tight oversampling: acceptance is p = 5/8 (eta=2) resp. 9/16 (eta=4).
  constexpr int BASE = ETA == 2 ? (8 * COUNT + 4) / 5 : (16 * COUNT + 8) / 9;
  constexpr int SAMPLE_COUNT = BASE + (int)isqrt_floor(2 * (int64_t)lambda * BASE) + 1;
  constexpr int BIT_COUNT = WIDTH * SAMPLE_COUNT;

  SharePair<nP, WIDTH, false> out;
  if constexpr (COUNT == 0)
    return out;

  const uint64_t and0 = bk.gmw().ands_evaluated;
  SharePair<nP, 1, true> dabits = unsigned_edabits<nP, 1>(bk, party, dealer, BIT_COUNT);

  // Gather the per-candidate bit shares.
  using Shares = emp::wrk::ShareVec<nP>;
  auto col = [&](int k) {
    Shares v((size_t)SAMPLE_COUNT);
    for (int j = 0; j < SAMPLE_COUNT; ++j)
      v[(size_t)j] = dabits.two_share[(size_t)WIDTH * j + k];
    return v;
  };
  auto Or = [&](const Shares& a, const Shares& b) {
    const auto ab = bk.gmw().and_batch(a, b);
    Shares o((size_t)SAMPLE_COUNT);
    for (int j = 0; j < SAMPLE_COUNT; ++j)
      o[(size_t)j] = emp::wrk::xor_share(emp::wrk::xor_share(a[(size_t)j], b[(size_t)j]), ab[(size_t)j]);
    return o;
  };

  Shares bad;
  if constexpr (ETA == 2) {
    // bad = b2 ^ (b1 & (b0 | b2))
    const Shares b0 = col(0), b1 = col(1), b2 = col(2);
    const Shares t = Or(b0, b2);
    const auto p = bk.gmw().and_batch(b1, t);
    bad.resize((size_t)SAMPLE_COUNT);
    for (int j = 0; j < SAMPLE_COUNT; ++j)
      bad[(size_t)j] = emp::wrk::xor_share(b2[(size_t)j], p[(size_t)j]);
  } else {
    // bad = b3 ^ (b2 & (b3 | b1 | b0))
    const Shares b0 = col(0), b1 = col(1), b2 = col(2), b3 = col(3);
    const Shares o = Or(Or(b0, b1), b3);
    const auto p = bk.gmw().and_batch(b2, o);
    bad.resize((size_t)SAMPLE_COUNT);
    for (int j = 0; j < SAMPLE_COUNT; ++j)
      bad[(size_t)j] = emp::wrk::xor_share(b3[(size_t)j], p[(size_t)j]);
  }

  const std::vector<uint8_t> opened_bad = bk.gmw().open(bad, "smallnorm-reject");
  std::vector<int> valid_index;
  valid_index.reserve((size_t)SAMPLE_COUNT);
  for (int j = 0; j < SAMPLE_COUNT; ++j)
    if (!opened_bad[(size_t)j])
      valid_index.push_back(j);
  if ((int)valid_index.size() < COUNT)
    std::fprintf(stderr, "rand DEBUG p%d: accepted %d of %d sampled, needed %d\n", party,
                 (int)valid_index.size(), SAMPLE_COUNT, COUNT);
  emp::expecting((int)valid_index.size() >= COUNT,
                 "rand: oversampled batch contained too few small-norm values");

  g_edabit_ands += bk.gmw().ands_evaluated - and0;

  // Keep the first COUNT accepted candidates; Boolean and arithmetic shares are
  // copied using the same source index, preserving each daBit-derived pair.
  out.fq_share.resize((size_t)COUNT);
  out.ring_share.resize((size_t)COUNT);
  out.two_share.resize((size_t)WIDTH * COUNT);
  for (int i = 0; i < COUNT; ++i) {
    const int src = valid_index[i];
    const int base = WIDTH * src;
    FqShare<nP> f = dabits.fq_share[base] + dabits.fq_share[base + 1] * 2;
    RingShare<nP> rs = dabits.ring_share[base] + dabits.ring_share[base + 1] * 2;
    if constexpr (ETA == 2) {
      f = f + dabits.fq_share[base + 2] * (uint32_t)(Q - 4);
      rs = rs - dabits.ring_share[base + 2] * 4;
    } else {
      f = f + dabits.fq_share[base + 2] * 4 + dabits.fq_share[base + 3] * (uint32_t)(Q - 8);
      rs = rs + dabits.ring_share[base + 2] * 4 - dabits.ring_share[base + 3] * 8;
    }
    out.fq_share[i] = f;
    out.ring_share[i] = rs;
    for (int k = 0; k < WIDTH; ++k)
      out.two_share[WIDTH * i + k] = dabits.two_share[base + k];
  }

  return out;
}

} // namespace mldsa
#endif // MLDSA_RAND_H
