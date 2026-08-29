// src/protocol/rand.h — F_smallnormpoly using daBits and emp-ag.
//
// For each output coefficient this returns matching authenticated shares:
//   q_share[i]                         = <r_i>_q over F_Q
//   two_share[BIT_WIDTH*i + k]         = <bit_k(r_i)>_2, little-endian
// where BIT_WIDTH is 3 for eta=2 and 4 for eta=4. The bit representation is
// two's complement, so the arithmetic value is
//   eta=2: b0 + 2*b1 - 4*b2
//   eta=4: b0 + 2*b1 + 4*b2 - 8*b3.
//
// The small-norm membership predicate is evaluated by an emp-ag circuit.
#ifndef MLDSA_RAND_H
#define MLDSA_RAND_H

#include "edabits.h"
#include <emp-ag/emp-ag.h>

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace mldsa {

inline AuthShare rand_scale(const AuthShare& x, uint32_t c) {
  return {fq_mul(x.val, c), fq_mul(x.mac, c)};
}

constexpr int64_t isqrt_floor(int64_t v) { // floor(sqrt(v)), constexpr-safe
  int64_t r = 0;
  while ((r + 1) * (r + 1) <= v)
    ++r;
  return r;
}

template <int nP, int ETA, int COUNT>
inline SharePair<nP, (ETA == 2 ? 3 : 4), false> rand_edabits(emp::AGMPCSession<nP>& sess, int party,
                                                             FakeDealer<nP>& dealer) {
  static_assert(ETA == 2 || ETA == 4, "rand: ETA must be 2 or 4");
  static_assert(COUNT >= 0, "rand: COUNT must be non-negative");

  constexpr int WIDTH = ETA == 2 ? 3 : 4;
  // Tight oversampling: acceptance is p = 5/8 (eta=2) resp. 9/16 (eta=4), so
  // BASE = COUNT/p draws suffice in expectation. A fixed additive +lambda is
  // NOT enough slack at this tight ratio (the binomial sd is ~sqrt(BASE)/2 >
  // lambda already at COUNT=1024); by Hoeffding a shortfall of t accepted
  // values happens w.p. <= exp(-2t^2/S), so t = sqrt(lambda*ln2*S/2) gives
  // 2^-lambda, and isqrt(2*lambda*BASE) covers t/p for both etas with margin.
  constexpr int BASE = ETA == 2 ? (8 * COUNT + 4) / 5 : (16 * COUNT + 8) / 9;
  constexpr int SAMPLE_COUNT = BASE + (int)isqrt_floor(2 * (int64_t)lambda * BASE) + 1;
  constexpr int BIT_COUNT = WIDTH * SAMPLE_COUNT;

  SharePair<nP, WIDTH, false> out;
  if constexpr (COUNT == 0)
    return out;

  SharePair<nP, 1, true> dabits = unsigned_edabits<nP, 1>(sess, party, dealer, BIT_COUNT);

  // Adopt every candidate bit in one AG input-label exchange, then evaluate
  // all rejection predicates in the same recorded circuit chunk.
  using Ctx = typename emp::AGMPCSession<nP>::ctx_t;
  using Bit = emp::Bit_T<Ctx>;
  using AllBits = emp::UInt_T<Ctx, BIT_COUNT>;
  AllBits bits = sess.template adopt_authenticated_input<AllBits>(dabits.two_share.data());
  std::vector<Bit> bad;
  bad.reserve((size_t)SAMPLE_COUNT);
  if constexpr (ETA == 2) {
    for (int j = 0; j < SAMPLE_COUNT; ++j) {
      const int base = 3 * j;
      bad.push_back(bits[base + 2] ^ (bits[base + 1] & (bits[base] | bits[base + 2])));
    }
  } else {
    for (int j = 0; j < SAMPLE_COUNT; ++j) {
      const int base = 4 * j;
      bad.push_back(bits[base + 3] ^
                    (bits[base + 2] & (bits[base + 3] | bits[base + 1] | bits[base])));
    }
  }

  std::vector<int> valid_index;
  valid_index.reserve((size_t)SAMPLE_COUNT);
  for (int j = 0; j < SAMPLE_COUNT; ++j) {
    const auto opened_bad = sess.reveal(bad[j], emp::PUBLIC);
    emp::expecting(opened_bad.has_value(), "rand: public predicate reveal failed");
    if (!opened_bad.value())
      valid_index.push_back(j);
  }
  if ((int)valid_index.size() < COUNT)
    std::fprintf(stderr, "rand DEBUG p%d: accepted %d of %d sampled, needed %d\n", party,
                 (int)valid_index.size(), SAMPLE_COUNT, COUNT);
  emp::expecting((int)valid_index.size() >= COUNT,
                 "rand: oversampled batch contained too few small-norm values");

  // Keep the first COUNT accepted candidates. Boolean and arithmetic shares
  // are copied using the same source index, preserving each daBit-derived pair.
  out.q_share.resize((size_t)COUNT);
  out.two_share.resize((size_t)WIDTH * COUNT);
  for (int i = 0; i < COUNT; ++i) {
    const int src = valid_index[i];
    const int base = WIDTH * src;
    AuthShare q = dabits.q_share[base];
    q = q + rand_scale(dabits.q_share[base + 1], 2);
    if constexpr (ETA == 2) {
      q = q + rand_scale(dabits.q_share[base + 2], Q - 4);
    } else {
      q = q + rand_scale(dabits.q_share[base + 2], 4);
      q = q + rand_scale(dabits.q_share[base + 3], Q - 8);
    }
    out.q_share[i] = q;
    for (int k = 0; k < WIDTH; ++k)
      out.two_share[WIDTH * i + k] = dabits.two_share[base + k];
  }

  return out;
}

} // namespace mldsa
#endif // MLDSA_RAND_H
