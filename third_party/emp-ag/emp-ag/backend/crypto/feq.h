#ifndef EMP_AG_BACKEND_CRYPTO_FEQ_H__
#define EMP_AG_BACKEND_CRYPTO_FEQ_H__
//
// F_eq (rush-safe equality test, eprint 2018/578) + the per-pair Fiat-Shamir coin
// machinery that derives challenges from the channel transcript. Network-stateful
// (takes NetIO / NetIOMP channels). Builds on the pure random-oracle hash fs_coin
// (defined at the top of this file).
//
#include "emp-ag/backend/netmp.h"
#include "emp-ag/backend/crypto/fs_hash.h"  // FsHash (FS/check digests)
#include "emp-tool/runtime/crypto/ro.h"   // RO (fs_coin random-oracle hash)
#include "emp-tool/emp-tool.h"

namespace emp::ag {

// Pure random-oracle helper: H(domain ‖ digest_a ‖ digest_b) over two channel-transcript
// digests. No transport, no protocol state. BYTE-IDENTITY CONTRACT: digest_a is absorbed
// BEFORE digest_b -- callers pass (ios[peer], ios2[peer]) in that order; swapping breaks
// the transcript. No broadcast, no extra rounds (it reuses the channel transcript the COT
// already hashes).
inline block fs_coin(const char *domain, block digest_a, block digest_b) {
  return RO(domain, zero_block).absorb(digest_a).absorb(digest_b).squeeze_block();
}

// Enable the IOChannel transcript hash on every pairwise channel (send_first =
// smaller-party convention). Guarded, so when the COT already enabled FS on a channel
// (SoftSpoken/Ferret do) this is a no-op and adds no hashing cost.
template <int nP>
inline void enable_fs_pairs(NetIOMP<nP> *io, int party) {
  for (int peer = 1; peer <= nP; ++peer)
    if (peer != party) {
      if (!io->ios[peer]->fs_enabled())  io->ios[peer]->enable_fs(party < peer);
      if (!io->ios2[peer]->fs_enabled()) io->ios2[peer]->enable_fs(party < peer);
    }
}
// Squeeze the per-pair coin: both parties of the pair fold the SAME two directed-channel
// digests (ios[peer] before ios2[peer]) in the same order, so they squeeze the identical
// coin -- a tweak consistent between them.
template <int nP>
inline block fs_coin_pair(NetIOMP<nP> *io, int peer, const char *domain) {
  return fs_coin(domain, io->ios[peer]->get_digest(),
                 io->ios2[peer]->get_digest());
}

// Fold a canonical set of lane transcripts into one local Fiat-Shamir coin.
// `digest_at(lane, direction)` returns direction 0 (`ios`) or 1 (`ios2`).
// Lane count/index framing prevents concatenation ambiguity, and the session id
// keeps otherwise identical multi-lane transcripts in separate executions from
// sharing a coin.
template <class DigestAt>
inline block fs_coin_lane_set(const char *domain, block session_id, int lanes,
                              DigestAt &&digest_at) {
  expecting(lanes > 0, "fs_coin_lane_set: lane count must be positive");
  RO ro(domain, session_id);
  ro.absorb((uint64_t)lanes);
  for (int lane = 0; lane < lanes; ++lane) {
    ro.absorb((uint64_t)lane);
    ro.absorb(digest_at(lane, 0));
    ro.absorb(digest_at(lane, 1));
  }
  return ro.squeeze_block();
}

// Multi-lane pair coin. Every byte transferred on a nonzero mesh is committed
// before the challenge becomes locally knowable. The one-lane case deliberately
// retains fs_coin_pair's historical output and wire-trace contract.
template <int nP>
inline block fs_coin_pair_group(NetIOMP<nP> *io, NetIOMPGroup<nP> *group,
                                int peer, const char *domain,
                                block session_id) {
  const int lanes = group == nullptr ? 1 : group->size();
  if (lanes == 1) return fs_coin_pair<nP>(io, peer, domain);
  return fs_coin_lane_set(
      domain, session_id, lanes, [io, group, peer](int lane, int direction) {
        NetIOMP<nP> *mesh = lane == 0 ? io : &group->mesh(lane);
        return direction == 0 ? mesh->ios[peer]->get_digest()
                              : mesh->ios2[peer]->get_digest();
      });
}

// Symmetric two-party F_eq, split so a caller may carry the commitment on an
// existing flight.  The commitment binds its role, session, and purpose: role
// binding prevents a rushing peer from reflecting the honest commitment and
// then copying its opening.
struct FeqSymmetric2PC {
  char digest[FsHash::DIGEST_SIZE];
  block nonce;
  char commitment[FsHash::DIGEST_SIZE];
  char peer_commitment[FsHash::DIGEST_SIZE];
};

inline void feq_symmetric_commitment_2pc(
    char out[FsHash::DIGEST_SIZE], block session_id, const char *purpose,
    int committer, const char digest[FsHash::DIGEST_SIZE], block nonce) {
  RO("emp-ag:F_eq:commit", session_id)
      .absorb(purpose)
      .absorb((uint64_t)committer)
      .absorb(digest, FsHash::DIGEST_SIZE)
      .absorb(nonce)
      .squeeze_digest(out);
}

inline FeqSymmetric2PC feq_prepare_symmetric_2pc(
    PRG *prg, block session_id, const char *purpose, int party,
    const char digest[FsHash::DIGEST_SIZE]) {
  FeqSymmetric2PC state{};
  memcpy(state.digest, digest, FsHash::DIGEST_SIZE);
  prg->random_block(&state.nonce, 1);
  feq_symmetric_commitment_2pc(state.commitment, session_id, purpose, party,
                               state.digest, state.nonce);
  return state;
}

// Complete the second (opening) flight after both commitments are fixed.
inline void feq_finish_symmetric_2pc(
    NetIO *send_io, NetIO *recv_io, block session_id, const char *purpose,
    int party, const FeqSymmetric2PC &state, const char *fail_msg) {
  send_io->send_data(state.digest, FsHash::DIGEST_SIZE);
  send_io->send_data(&state.nonce, sizeof(block));
  send_io->flush();

  char peer_digest[FsHash::DIGEST_SIZE];
  block peer_nonce;
  recv_io->recv_data(peer_digest, FsHash::DIGEST_SIZE);
  recv_io->recv_data(&peer_nonce, sizeof(block));

  char expected[FsHash::DIGEST_SIZE];
  feq_symmetric_commitment_2pc(expected, session_id, purpose, 3 - party,
                               peer_digest, peer_nonce);
  expecting(memcmp(expected, state.peer_commitment, FsHash::DIGEST_SIZE) == 0,
            "F_eq: symmetric commit-open mismatch");
  expecting(memcmp(state.digest, peer_digest, FsHash::DIGEST_SIZE) == 0,
            fail_msg);
}

}  // namespace emp::ag
#endif  // EMP_AG_BACKEND_CRYPTO_FEQ_H__
