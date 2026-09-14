#ifndef EMP_AG_BACKEND_CRYPTO_FZERO_H__
#define EMP_AG_BACKEND_CRYPTO_FZERO_H__
//
// Π_FZero: shared-zero block-vector via pairwise seed-and-expand. Network-stateful
// (NetIOMP mesh + ThreadPool + PRG).
//
#include "emp-ag/backend/runtime.h"   // parallel_for, joinNclean
#include "emp-ag/backend/netmp.h"
#include "emp-tool/emp-tool.h"

namespace emp::ag {

using std::future;
using std::vector;

// Π_FZero: shared-zero block-vector via pairwise seed-and-expand under
// the ICM assumption on the AES-based PRG. Each peer-pair shares a
// λ-bit seed (the larger-indexed party samples and sends to the
// smaller); both expand the seed via PRG into n blocks and XOR into
// `out`. Each seed contributes to exactly two parties so Σ_p out^p[k]
// = 0 for every k. Communication: λ bits per pair (vs nP·λ for the
// straightforward zero-share). Caller is responsible for zero-init —
// the contributions are XORed in so this composes with caller buffers
// that already have content (e.g. TriplePool loads z directly into phi).
template <int nP>
void fzero_xor(NetIOMP<nP> *io, PRG *prg, ThreadPool *pool, int party,
               block *out, int n, BlockVec &scratch) {
  block seed_by_peer[nP + 1];
  prg->random_block(&seed_by_peer[1], nP);

  vector<future<void>> res;
  for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
    res.push_back(pool->enqueue([io, &seed_by_peer, peer, party]() {
      if (peer < party) {
        io->send_data(peer, &seed_by_peer[peer], sizeof(block));
        io->flush(peer);
      } else {
        io->recv_data(peer, &seed_by_peer[peer], sizeof(block));
      }
    }));
  }
  joinNclean(res);

  if ((int)scratch.size() < n) scratch.resize(n);   // caller-owned grow-only (no churn)
  // Expand + fold, range-split: CTR position c of a seed's stream is the same
  // block no matter which worker produces it (PRG::fork_at), so each range forks
  // every peer's expansion at its own offset — bit-identical to the sequential
  // per-peer expansion, with disjoint scratch/out slices per range. This was a
  // serial (n-1) x n-block expansion on the leaky-AND critical path.
  PRG base[nP + 1];
  for (int peer = 1; peer <= nP; ++peer) if (peer != party)
    base[peer].reseed(&seed_by_peer[peer]);
  block *sc = scratch.data();
  ag::parallel_for(pool, n, [&](int lo, int hi) {
    for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
      PRG p = base[peer].fork_at((uint64_t)lo);
      p.random_block(sc + lo, hi - lo);
      xorBlocks_arr(out + lo, out + lo, sc + lo, hi - lo);
    }
  });
}

}  // namespace emp::ag
#endif  // EMP_AG_BACKEND_CRYPTO_FZERO_H__
