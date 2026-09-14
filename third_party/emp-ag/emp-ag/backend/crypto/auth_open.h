#pragma once

// Local addition: stand-alone authenticated Boolean openings. Upstream
// preprocessing defers error detection to its garbling checker; this adapter
// verifies every claimed bit before using it in an authenticated GMW result.
#include "emp-ag/backend/crypto/fs_hash.h"
#include "emp-ag/backend/netmp.h"
#include "emp-ag/backend/runtime.h"
#include <array>
#include <climits>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace emp::ag {

inline void auth_open_hash_u64(FsHash &h, uint64_t value) {
  unsigned char bytes[8];
  for (int i = 0; i < 8; ++i) bytes[i] = static_cast<unsigned char>(value >> (8 * i));
  h.put(bytes, 8);
}
inline void auth_open_prefix(FsHash &h, block session, const std::string &domain,
                             uint64_t epoch, uint64_t count, int holder, int verifier) {
  static constexpr char name[] = "ThresholdMLDSA/authenticated-Boolean-open/v1";
  h.put(name, sizeof(name)); h.put(&session, sizeof(session));
  auth_open_hash_u64(h, domain.size()); h.put(domain.data(), domain.size());
  auth_open_hash_u64(h, epoch); auth_open_hash_u64(h, count);
  auth_open_hash_u64(h, holder); auth_open_hash_u64(h, verifier);
}

// views[p][j] is party p's declared XOR share. MAC(peer,j) is this holder's
// tag under peer's key; KEY(peer,j) authenticates that peer under local Delta.
// Send a hash of the entire MAC vector, not a reusable linear MAC checksum.
// Every honest party broadcasts its local check result and the same canonical
// opening view before anyone returns. A peer disconnect is a protocol abort.
template<int N, class Mac, class Key>
void authenticate_open_claims(NetIOMP<N> *io, ThreadPool *pool, int party,
    block Delta, block session, const std::string &domain, uint64_t epoch,
    size_t count, const unsigned char *const views[N + 1], Mac mac, Key key) {
  using Digest = std::array<unsigned char, FsHash::DIGEST_SIZE>;
  std::array<Digest, N + 1> sent{}, received{}, expected{};
  bool good = true;
  for (int peer = 1; peer <= N; ++peer) if (peer != party) {
    FsHash hs, he;
    auth_open_prefix(hs, session, domain, epoch, count, party, peer);
    auth_open_prefix(he, session, domain, epoch, count, peer, party);
    for (size_t j = 0; j < count; ++j) {
      const block m = mac(peer, j);
      good = good && views[peer][j] <= 1 && views[party][j] <= 1;
      const block e = key(peer, j) ^ (select_mask[views[peer][j] & 1] & Delta);
      hs.put(&m, sizeof(m)); he.put(&e, sizeof(e));
    }
    hs.digest(reinterpret_cast<char *>(sent[peer].data()));
    he.digest(reinterpret_cast<char *>(expected[peer].data()));
  }
  std::vector<std::future<void>> jobs;
  for (int peer = 1; peer <= N; ++peer) if (peer != party) {
    jobs.push_back(pool->enqueue([&, peer] {
      io->send_data(peer, sent[peer].data(), sent[peer].size()); io->flush(peer);
    }));
    jobs.push_back(pool->enqueue([&, peer] {
      io->recv_data(peer, received[peer].data(), received[peer].size());
    }));
  }
  joinNclean(jobs);
  for (int peer = 1; peer <= N; ++peer) if (peer != party)
    good = good && received[peer] == expected[peer];

  FsHash vh;
  auth_open_prefix(vh, session, domain, epoch, count, 0, 0);
  for (int p = 1; p <= N; ++p) vh.put(views[p], count);
  std::array<unsigned char, 1 + FsHash::DIGEST_SIZE> decision{};
  decision[0] = good ? 1 : 0;
  vh.digest(reinterpret_cast<char *>(decision.data() + 1));
  std::array<decltype(decision), N + 1> decisions{};
  for (int peer = 1; peer <= N; ++peer) if (peer != party) {
    jobs.push_back(pool->enqueue([&, peer] {
      io->send_data(peer, decision.data(), decision.size()); io->flush(peer);
    }));
    jobs.push_back(pool->enqueue([&, peer] {
      io->recv_data(peer, decisions[peer].data(), decisions[peer].size());
    }));
  }
  joinNclean(jobs);
  for (int peer = 1; peer <= N; ++peer) if (peer != party)
    good = good && decisions[peer][0] == 1 &&
        std::memcmp(decision.data() + 1, decisions[peer].data() + 1,
                    FsHash::DIGEST_SIZE) == 0;
  if (!good) throw std::runtime_error("authenticated Boolean open failed: " + domain);
}

template<int N, class Mac, class Key>
std::vector<uint8_t> checked_open_bits(NetIOMP<N> *io, ThreadPool *pool, int party,
    block Delta, block session, const std::string &domain, uint64_t epoch,
    const std::vector<uint8_t> &mine, Mac mac, Key key) {
  const size_t count = mine.size();
  if (count > static_cast<size_t>(INT_MAX)) throw std::runtime_error("Boolean open too large");
  std::array<std::vector<unsigned char>, N + 1> views;
  for (int p = 1; p <= N; ++p) views[p].resize(count);
  views[party].assign(mine.begin(), mine.end());
  std::vector<std::future<void>> jobs;
  for (int peer = 1; peer <= N; ++peer) if (peer != party) {
    jobs.push_back(pool->enqueue([&, peer] {
      io->send_bool(peer, reinterpret_cast<const bool *>(mine.data()), count); io->flush(peer);
    }));
    jobs.push_back(pool->enqueue([&, peer] {
      io->recv_bool(peer, reinterpret_cast<bool *>(views[peer].data()), count);
    }));
  }
  joinNclean(jobs);
  const unsigned char *ptrs[N + 1]{};
  for (int p = 1; p <= N; ++p) ptrs[p] = views[p].data();
  authenticate_open_claims(io, pool, party, Delta, session, domain, epoch, count,
                           ptrs, mac, key);
  std::vector<uint8_t> result(count, 0);
  for (int p = 1; p <= N; ++p)
    for (size_t j = 0; j < count; ++j) result[j] ^= views[p][j];
  return result;
}
}  // namespace emp::ag
