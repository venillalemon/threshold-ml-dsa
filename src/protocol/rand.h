// src/protocol/rand.h — F_smallnormpoly, sampled and shared directly by the dealer.
//
// Each coefficient is sampled uniformly in [-eta, eta] (no rejection sampling)
// and dealt as matching F_q and ring shares. s, e (KeyGen) and e_w (PrepSign)
// only ever use these arithmetic shares, so no Boolean half is produced and no
// GMW runs here -- these are correlated setup values the paper assumes.
#ifndef MLDSA_RAND_H
#define MLDSA_RAND_H

#include "edabits.h"

#include <cstdint>
#include <vector>

namespace mldsa {

// Retained for API compatibility with cost accounting; edaBits no longer spend
// GMW AND gates, so this stays zero.
inline uint64_t g_edabit_ands = 0;

template <int nP, int ETA, int COUNT>
inline SharePair<nP, (ETA == 2 ? 3 : 4), false> rand_edabits(Backend<nP>& /*bk*/, int /*party*/,
                                                             FakeDealer<nP>& dealer) {
  static_assert(ETA == 2 || ETA == 4, "rand: ETA must be 2 or 4");
  static_assert(COUNT >= 0, "rand: COUNT must be non-negative");
  constexpr int WIDTH = ETA == 2 ? 3 : 4;

  SharePair<nP, WIDTH, false> out;
  out.fq_share.resize((size_t)COUNT);
  out.ring_share.resize((size_t)COUNT);
  for (int i = 0; i < COUNT; ++i) {
    const int64_t val = (int64_t)dealer.draw_uniform(2 * (uint64_t)ETA + 1) - ETA; // [-eta, eta]
    const uint32_t vq = (uint32_t)(((val % Q) + Q) % Q);
    out.fq_share[(size_t)i] = dealer.deal_fq(vq);
    out.ring_share[(size_t)i] = dealer.deal_ring(ring_from_i64(val));
  }
  return out;
}

} // namespace mldsa
#endif // MLDSA_RAND_H
