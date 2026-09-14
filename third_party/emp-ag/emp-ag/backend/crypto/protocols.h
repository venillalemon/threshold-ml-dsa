#ifndef EMP_AG_BACKEND_CRYPTO_PROTOCOLS_H__
#define EMP_AG_BACKEND_CRYPTO_PROTOCOLS_H__
//
// Subprotocol primitives over the NetIOMP mesh: the equivocation-safe echo
// broadcast (EchoBC), the multiparty commit-then-open coin (sampleRandom), the
// commit-open shared-zero test (check_shared_zero), and the DEBUG-ONLY
// check_MAC. Network-stateful but protocol-agnostic: nothing here knows about
// garbling, triples, or the engine.
//
#include "emp-ag/backend/netmp.h"
#include "emp-ag/backend/crypto/fs_hash.h"  // FsHash (FS/check digests)
#include "emp-ag/backend/runtime.h"
#include "emp-tool/runtime/crypto/ro.h"
#include "emp-tool/emp-tool.h"
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace emp::ag {

// Commitment shared by the preprocessing commit/open checks. Binding the
// session, purpose, and committer prevents a rushing corrupt party from copying
// an honest party's commitment and later reflecting its opening under another
// party index. Payload length is framed explicitly before the bytes.
inline void role_bound_commitment(
    char out[FsHash::DIGEST_SIZE], block session_id, const char *purpose,
    int committer, const void *payload, size_t payload_len) {
  RO("emp-ag:protocol-commit", session_id)
      .absorb(purpose)
      .absorb((uint64_t)committer)
      .absorb((uint64_t)payload_len)
      .absorb(payload, payload_len)
      .squeeze_digest(out);
}

// Echo broadcast with deferred finalize. Each all_bcast does an insecure
// pairwise send (one round) and folds the received view into a rolling
// transcript hash in canonical order. finalize() exchanges transcript
// digests pairwise and aborts on any mismatch — one round suffices because
// honest parties compare received digests against their own. Any logical
// broadcast routed through all_bcast is protected against equivocation.
template <int nP>
class EchoBC {
public:
  FsHash h;
  NetIOMP<nP> *io;
  ThreadPool *pool;
  int party;

  EchoBC(NetIOMP<nP> *io_in, ThreadPool *pool_in, int party_in)
      : io(io_in), pool(pool_in), party(party_in) {}

  // Broadcast `len` elements of type T per party. view[p] must point to a
  // buffer of at least `len` T's for every p; view[party] is populated from
  // my_v, and view[p] for p != party is filled from the wire. The full
  // view (nP * len * sizeof(T) bytes) is folded into the transcript hash
  // in party-index order.
  template <typename T>
  void all_bcast(const T *my_v, int len, T *view[nP + 1]) {
    memcpy(view[party], my_v, sizeof(T) * len);
    int bytes = (int)(sizeof(T) * len);
    vector<future<void>> res;
    for (int i = 1; i <= nP; ++i) for (int j = 1; j <= nP; ++j)
      if ((i < j) and (i == party or j == party)) {
        int p2 = i + j - party;
        res.push_back(pool->enqueue([this, my_v, bytes, p2]() {
          io->send_data(p2, my_v, bytes);
          io->flush(p2);
        }));
        res.push_back(pool->enqueue([this, &view, bytes, p2]() {
          io->recv_data(p2, view[p2], bytes);
        }));
      }
    joinNclean(res);
    for (int p = 1; p <= nP; ++p) h.put(view[p], bytes);
  }

  // Single-element convenience overload: view is a flat [nP+1] stack array.
  template <typename T>
  void all_bcast(const T &my_v, T view[nP + 1]) {
    T *ptrs[nP + 1];
    for (int p = 1; p <= nP; ++p) ptrs[p] = &view[p];
    all_bcast(&my_v, 1, ptrs);
  }

  void finalize() {
    char d_me[FsHash::DIGEST_SIZE];
    h.digest(d_me);
    char recv[nP + 1][FsHash::DIGEST_SIZE];
    vector<future<void>> res;
    for (int i = 1; i <= nP; ++i) for (int j = 1; j <= nP; ++j)
      if ((i < j) and (i == party or j == party)) {
        int p2 = i + j - party;
        res.push_back(pool->enqueue([this, &d_me, p2]() {
          io->send_data(p2, d_me, FsHash::DIGEST_SIZE);
          io->flush(p2);
        }));
        res.push_back(pool->enqueue([this, &recv, p2]() {
          io->recv_data(p2, recv[p2], FsHash::DIGEST_SIZE);
        }));
      }
    joinNclean(res);
    for (int p = 1; p <= nP; ++p) if (p != party)
      expecting(memcmp(d_me, recv[p], FsHash::DIGEST_SIZE) == 0,
                "EchoBC finalize: transcript divergence\n");
  }
};

template <int nP>
block sampleRandom(NetIOMP<nP> *io, PRG *prg, ThreadPool *pool, int party,
                   block session_id = zero_block) {
  vector<future<void>> res;
  vector<future<bool>> res2;
  char(*dgst)[FsHash::DIGEST_SIZE] = new char[nP + 1][FsHash::DIGEST_SIZE];
  char echo[nP + 1][FsHash::DIGEST_SIZE];
  block *S = new block[nP + 1];
  prg->random_block(&S[party], 1);
  auto commit = [session_id](char *out, int sender, const block &seed) {
    RO("emp-ag:sample-random:commit", session_id)
        .absorb((uint64_t)sender)
        .absorb(seed)
        .squeeze_digest(out);
  };
  commit(dgst[party], party, S[party]);

  for (int i = 1; i <= nP; ++i)
    for (int j = 1; j <= nP; ++j)
      if ((i < j) and (i == party or j == party)) {
        int party2 = i + j - party;
        res.push_back(pool->enqueue([dgst, io, party, party2]() {
          io->send_data(party2, dgst[party], FsHash::DIGEST_SIZE);
          io->flush(party2);
          io->recv_data(party2, dgst[party2], FsHash::DIGEST_SIZE);
        }));
      }
  joinNclean(res);

  // Echo C1 without another round: carry the round-1 view digest alongside S.
  FsHash::hash_once(echo[party], dgst[1], nP * FsHash::DIGEST_SIZE);
  for (int i = 1; i <= nP; ++i)
    for (int j = 1; j <= nP; ++j)
      if ((i < j) and (i == party or j == party)) {
        int party2 = i + j - party;
        res2.push_back(pool->enqueue([io, S, dgst, &echo, commit, party, party2]() -> bool {
          io->send_data(party2, &S[party], sizeof(block));
          io->send_data(party2, echo[party], FsHash::DIGEST_SIZE);
          io->flush(party2);
          io->recv_data(party2, &S[party2], sizeof(block));
          io->recv_data(party2, echo[party2], FsHash::DIGEST_SIZE);
          char tmp[FsHash::DIGEST_SIZE];
          commit(tmp, party2, S[party2]);
          return memcmp(tmp, dgst[party2], FsHash::DIGEST_SIZE) != 0 ||
                 memcmp(echo[party], echo[party2], FsHash::DIGEST_SIZE) != 0;
        }));
      }
  bool cheat = joinNcleanCheat(res2);
  expecting(!cheat, "cheat in sampleRandom");
  for (int i = 2; i <= nP; ++i) S[1] = S[1] ^ S[i];
  block result = S[1];
  delete[] S;
  delete[] dgst;
  return result;
}

// Commit-open shared-zero test over the mesh (nP > 2): each party holds a block
// m_i; verify XOR_i m_i = 0 without revealing anything else. Rush-safe
// (commit-before-open: each commitment binds session, purpose, committer,
// m, and r; commitments are echo-broadcast first and openings second) and
// equivocation-safe (both rounds
// ride one EchoBC; finalize() confirms every party saw the same view). Aborts
// via expecting(..., fail_msg) on a nonzero sum and
// expecting(..., commit_fail_msg) on a
// commit-open mismatch — both caller-supplied, so the abort text names the
// caller's check, not this primitive.
template <int nP>
void check_shared_zero(NetIOMP<nP> *io, ThreadPool *pool, PRG *prg, int party,
                       const block &m_i, block session_id, const char *purpose,
                       const char *fail_msg,
                       const char *commit_fail_msg) {
  struct OpenT { block m, r; };
  struct DigestT { char d[FsHash::DIGEST_SIZE]; };
  OpenT mine;
  mine.m = m_i;
  prg->random_block(&mine.r, 1);
  DigestT com_mine;
  role_bound_commitment(com_mine.d, session_id, purpose, party,
                        &mine, sizeof(OpenT));              // r keeps it hiding

  EchoBC<nP> echo(io, pool, party);
  DigestT com_view[nP + 1];
  echo.all_bcast(com_mine, com_view);                    // round 1: commitments
  OpenT open_view[nP + 1];
  echo.all_bcast(mine, open_view);                       // round 2: openings

  block sum = zero_block;
  for (int p = 1; p <= nP; ++p) {
    DigestT chk;
    role_bound_commitment(chk.d, session_id, purpose, p,
                          &open_view[p], sizeof(OpenT));
    expecting(memcmp(chk.d, com_view[p].d, FsHash::DIGEST_SIZE) == 0,
              commit_fail_msg);
    sum = sum ^ open_view[p].m;
  }
  expecting(cmpBlock(&sum, &zero_block, 1), fail_msg);
  echo.finalize();                                       // transcript agreement
}

// DEBUG-ONLY MAC-consistency check. WARNING: broadcasts each party's Delta in the CLEAR (line
// below) -- a malicious party could then forge any MAC and defeat the c_gamma MACCheck. The
// DEFINITIONS below are compiled ONLY under EMP_AG_DEBUG_CHECKS, so a production build cannot
// even instantiate them; the white-box tests (auth_share, three_party_iknp) opt in. Must NEVER
// ship enabled.
#ifdef EMP_AG_DEBUG_CHECKS
template <int nP>
void check_MAC(NetIOMP<nP> *io, block *MAC[nP + 1], block *KEY[nP + 1], bool *r,
               block Delta, int length, int party) {
  block *tmp = new block[length];
  block tD;
  for (int i = 1; i <= nP; ++i)
    for (int j = 1; j <= nP; ++j)
      if (i < j) {
        if (party == i) {
          io->send_data(j, &Delta, sizeof(block));
          io->send_data(j, KEY[j], sizeof(block) * length);
          io->flush(j);
        } else if (party == j) {
          io->recv_data(i, &tD, sizeof(block));
          io->recv_data(i, tmp, sizeof(block) * length);
          for (int k = 0; k < length; ++k) {
            if (r[k])
              tmp[k] = tmp[k] ^ tD;
          }
          expecting(cmpBlock(MAC[i], tmp, length), "check_MAC failed!");
        }
      }
  delete[] tmp;
  if (party == 1)
    cerr << "check_MAC pass!\n" << flush;
}

template <int nP>
void check_MAC(NetIOMP<nP> *io, BlockVec MAC[nP + 1],
               BlockVec KEY[nP + 1], std::vector<unsigned char> &r, block Delta,
               int length, int party) {
  block *MAC_p[nP + 1], *KEY_p[nP + 1];
  for (int i = 1; i <= nP; ++i) if (i != party) {
    MAC_p[i] = MAC[i].data();
    KEY_p[i] = KEY[i].data();
  }
  check_MAC<nP>(io, MAC_p, KEY_p, (bool *)r.data(), Delta, length, party);
}

// r-less overload: derive share bits from bit0(MAC[any peer]). Valid when
// the pool invariant holds (bit0(K)=0, bit0(Δ)=1 ⇒ bit0(M) = x consistently
// across peers).
template <int nP>
void check_MAC(NetIOMP<nP> *io, BlockVec MAC[nP + 1],
               BlockVec KEY[nP + 1], block Delta, int length, int party) {
  int any_peer = (party == 1) ? 2 : 1;
  std::vector<unsigned char> r(length);
  for (int k = 0; k < length; ++k)
    r[k] = (unsigned char)getLSB(MAC[any_peer][k]);
  check_MAC<nP>(io, MAC, KEY, r, Delta, length, party);
}
#endif  // EMP_AG_DEBUG_CHECKS
}  // namespace emp::ag
#endif // EMP_AG_BACKEND_CRYPTO_PROTOCOLS_H__
