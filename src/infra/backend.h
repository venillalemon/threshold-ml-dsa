// src/infra/backend.h — the mixed-security boolean backend session.
//
// This replaces the single KRRW AGMPCSession of the earlier prototype with the
// split the paper actually asks for:
//   * offline boolean circuits (C_pre / Phase 1) run under authenticated GMW
//     (emp::gmw::Gmw, real malicious-COT TinyOT-style triples);
//   * the online circuit (C_post / Phase 2) is a WRK17b four-row garbled circuit
//     (emp::wrk), garbled and delivered — with its output masks opened to P1 —
//     entirely in the offline phase, so signing costs exactly two post-challenge
//     flights (three in slot mode, for the routed secret input).
//
// One NetIOMP<nP> is shared by three consumers, used strictly sequentially:
//   * the GMW engine (its own TriplePool owns the COT setup);
//   * the WRK offline delivery (rows + fixed labels + output masks to P1);
//   * the arithmetic BDOZ opens of spdz.h (open_fq_checked / open_ring_* ).
// The arithmetic layer (FakeDealer) keeps its own slopes; the boolean layer
// lives entirely in emp::wrk::AuthShare under this session's single GMW Delta.
#ifndef MLDSA_BACKEND_H
#define MLDSA_BACKEND_H

#include "spdz.h" // mldsa::deterministic_delta (shared with FakeDealer)

#include <emp-ag/gmw.h>
#include <emp-ag/wrk.h>
#include <emp-ag/backend/netmp.h>

#include <cstdint>
#include <cstring>
#include <memory>

namespace mldsa {

// The backend's Delta is the SAME deterministic per-party key the FakeDealer
// reconstructs (spdz.h), so dealer-fabricated Boolean edaBit shares authenticate
// under it. Demo cheat: a real deployment keeps Delta private.

// A fixed, party-independent session id. Every party derives the identical
// block with no communication. This is a demo stand-in for a fresh unpredictable
// per-execution id (like the FakeDealer's shared seed); a real deployment binds
// it to the agreed transcript context.
inline emp::block demo_session_id(uint64_t seed = 0x7468726573686f6cULL /* "threshol" */) {
  return emp::makeBlock(seed, 0xd6d6c6473610000ULL);
}

// The boolean backend. Owns the mesh, the thread pool reference, and the GMW
// engine; hands out io() for the arithmetic BDOZ opens and gmw() for offline
// circuits. WRK online state is built by the wrk_phase2.h glue against gmw().
template <int nP> class Backend {
public:
  using Share = emp::wrk::AuthShare<nP>;
  using Shares = emp::wrk::ShareVec<nP>;

  Backend(int party, int port, ThreadPool* pool, int triple_ssp = 80)
      : party_(party), pool_(pool), io_(party, port),
        delta_(deterministic_delta<nP>(party)), sid_(demo_session_id()),
        gmw_(&io_, pool, party, sid_, delta_, triple_ssp) {}

  Backend(const Backend&) = delete;
  Backend& operator=(const Backend&) = delete;

  int party() const { return party_; }
  emp::block delta() const { return delta_; }
  emp::block session_id() const { return sid_; }
  emp::NetIOMP<nP>& io() { return io_; }
  ThreadPool* pool() { return pool_; }
  emp::gmw::Gmw<nP>& gmw() { return gmw_; }

  // Total bytes this party has moved on the mesh (sent + received).
  int64_t comm_bytes() { return io_.count(); }

private:
  int party_;
  ThreadPool* pool_;
  emp::NetIOMP<nP> io_;
  emp::block delta_;
  emp::block sid_;
  emp::gmw::Gmw<nP> gmw_;
};

} // namespace mldsa

#endif // MLDSA_BACKEND_H
