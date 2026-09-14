#ifndef EMP_AG_BACKEND_RUNTIME_H__
#define EMP_AG_BACKEND_RUNTIME_H__
//
// Execution-infra shared by the backend — pool sizing, the lane partition (the
// AGREED wire-shaping constants kLaneCRef / lane_range), width-splits, and the
// SPSC chunk pipe. Pure scheduling/partitioning state: no network mesh, no
// protocol state.
//
#include "emp-tool/emp-tool.h"
#include <condition_variable>
#include <future>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace emp::ag {

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::future;
using std::max;
using std::vector;

// The LaAND row window carries ceil((ssp+2)/8) bytes in one 128-bit block.
// Hence ssp=126 is the largest value whose advertised number of checked bits
// fits the construction; non-positive values are not meaningful security
// parameters. Keep the pure predicate available to configuration tests.
inline constexpr int kMinStatisticalSecurity = 1;
inline constexpr int kMaxStatisticalSecurity = 126;
inline constexpr bool valid_statistical_security(int ssp) {
  return ssp >= kMinStatisticalSecurity && ssp <= kMaxStatisticalSecurity;
}

inline int checked_statistical_security(int ssp) {
  expecting(valid_statistical_security(ssp), "AGMPC: ssp must be in [1, 126]");
  return ssp;
}

// Interactive protocol phases enqueue mutually-dependent send/receive work.
// A null/empty pool is immediately invalid, and a single worker can block in
// one direction while the matching progress task remains queued.
inline void validate_protocol_runtime(ThreadPool *pool, int party, int nP) {
  expecting(party >= 1 && party <= nP, "AGMPC: party must be in [1, nP]");
  expecting(pool != nullptr, "AGMPC: ThreadPool must not be null");
  expecting(pool->size() >= 2, "AGMPC: ThreadPool must have at least two workers");
}

// Single-producer / single-consumer bounded ring of pre-allocated chunk buffers.
// The producer (a garble pass) fills slots and publishes them; the consumer (a
// drain thread) ships each slot and releases it, so per-chunk compute overlaps
// per-chunk network I/O. producer_slot() blocks while the ring is full;
// consumer_slot() blocks while it is empty and returns nullptr once the producer
// has closed and the ring is drained. Slots are constructed once via the Init
// callback so a chunk's buffers are sized to the per-call cap up front (no
// allocation in the hot path).
template <typename T>
class chunk_pipe {
public:
  template <typename Init>
  chunk_pipe(size_t depth, Init init) : slots_(depth) {
    for (auto &s : slots_) init(s);
  }

  T &producer_slot() {
    std::unique_lock<std::mutex> lk(mu_);
    cv_free_.wait(lk, [&] { return count_ < slots_.size(); });
    return slots_[(consumer_head_ + count_) % slots_.size()];
  }
  void producer_publish() {
    {
      std::lock_guard<std::mutex> lk(mu_);
      ++count_;
    }
    cv_full_.notify_one();
  }
  void producer_close() {
    {
      std::lock_guard<std::mutex> lk(mu_);
      closed_ = true;
    }
    cv_full_.notify_one();
  }

  T *consumer_slot() {
    std::unique_lock<std::mutex> lk(mu_);
    cv_full_.wait(lk, [&] { return count_ > 0 || closed_; });
    if (count_ == 0) return nullptr;            // drained + closed
    return &slots_[consumer_head_];
  }
  void consumer_release() {
    {
      std::lock_guard<std::mutex> lk(mu_);
      consumer_head_ = (consumer_head_ + 1) % slots_.size();
      --count_;
    }
    cv_free_.notify_one();
  }

private:
  std::vector<T> slots_;
  size_t consumer_head_ = 0;
  size_t count_ = 0;
  bool closed_ = false;
  std::mutex mu_;
  std::condition_variable cv_full_, cv_free_;
};

// Default ThreadPool size: the LOCAL core count (with a two-worker progress
// floor), independent of nP.
// The party count sets how many work items are
// ENQUEUED (2*(nP-1) COT streams, x L lanes), not the pool size: at small nP the
// lanes fill the cores; at large nP the surplus stream tasks run in FIFO waves
// over the core-sized pool (the NIC saturates well below the core count, so waves
// are throughput-neutral). Compute splits (parallel_for) self-cap at pool.size()
// == cores, so they never over-split. Pool size is a LOCAL resource and may
// differ across parties; the wire-shaping L(nP) is agreed separately.
//   (Supersedes the min(2*(nP-1), cores) cap of commit 1007058, which left cores
//   idle at small nP and mis-capped the I/O concurrency at large nP.)
inline int default_pool_size([[maybe_unused]] int nP) {
  int hw = (int)std::thread::hardware_concurrency();
  return std::max(2, hw > 0 ? hw : 2 * (nP - 1));
}

// Production lane count L(nP) for SoftSpoken lane-parallel COT.
// Each party runs
// 2*(nP-1) base COT streams; L lanes multiply that to fill the core budget with
// COT-stream concurrency (the binding constraint on the draw+seed phase --
// pool size alone cannot parallelize it, see the thread-sweep finding).
//
// kLaneCRef is the TARGET total lane-stream concurrency 2*(nP-1)*L; it is an
// AGREED value that shapes the wire (both parties must pin it identically, like
// SoftSpoken's k) -- NOT derived from local cores, which differ across machines.
// The default C_ref=8 targets the small-nP knee: it gives L(2)=4 and L(3)=2,
// which sit at (or within noise of) the per-nP COT-stream saturation point --
// a sharper nP=2 knee than C_ref=16 would.
inline constexpr int kLaneCRef = 8;   // agreed target lane-stream concurrency
inline constexpr int kLaneMax  = 4;   // socket/memory cap on L
inline constexpr int lanes_for(int nP) {
  int L = kLaneCRef / (2 * (nP - 1));
  if (L < 1) L = 1;
  if (L > kLaneMax) L = kLaneMax;
  return L;   // L(2)=4, L(3)=2, L(4)=1, ...
}

// Contiguous slice [lo, hi) of [0, n) for lane t of T lanes: an even split with the
// first (n mod T) lanes taking one extra. THE single source of the lane partition —
// COT (AuthSharePool) and every triple_pool bulk-open fan-out slice through this, so
// both ends compute the same boundaries (R1). Pure function of (t, n, T).
inline std::pair<int, int> lane_range(int t, int n, int T) {
  const int base = n / T, rem = n % T;
  const int lo = t * base + std::min(t, rem);
  return {lo, lo + base + (t < rem ? 1 : 0)};
}

// Peer visitation order for the pairwise all-to-all fan-outs (COT etc.). Default is
// plain index order (1..nP skip self) — byte-identical to the historical fan-out. Under
// EMP_AG_PEER_XOR, and when nP is a power of 2, the XOR 1-FACTORIZATION: round r=1..nP-1
// pairs `party` with `(party-1)^r+1`. Its point is that BOTH ends of pair (i,j) land in
// the SAME round r=(i-1)^(j-1), so they are scheduled simultaneously. In the WAVE regime
// (pool < 2(nP-1) tasks, i.e. large nP or few cores) that keeps every party busy each
// round — a perfect matching per round (deadlock-trivial) — instead of the peer-index
// WAVEFRONT where a high party stalls waiting for a low peer whose turn has not come.
// For nP=2^d the powers-of-2 r are the hypercube dimensions; the rest fill the complete
// graph. Reordering is BYTE-IDENTICAL: each pair's COT is per-channel deterministic, so
// only the schedule changes, not the wire. Non-power-of-2 nP falls back to index order
// (circle-method 1-factorization is the general TODO).
inline std::vector<int> peer_order(int party, int nP) {
  std::vector<int> peers;
  peers.reserve(nP > 0 ? (size_t)(nP - 1) : 0);
#ifdef EMP_AG_PEER_XOR
  if (nP >= 2 && (nP & (nP - 1)) == 0) {   // power of 2
    for (int r = 1; r < nP; ++r) peers.push_back(((party - 1) ^ r) + 1);
    return peers;
  }
#endif
  for (int p = 1; p <= nP; ++p) if (p != party) peers.push_back(p);
  return peers;
}

// Width-split [0, n) across the pool: body(a, b) processes the half-open range
// [a, b). Caller guarantees the ranges' writes are disjoint. One join. Used for
// the serial single-thread loops on the protocol's critical path (Beaver fold,
// x/y reveal math, the RLC and coeff expansions, the d-combine) -- pure
// critical-path time on otherwise-idle budgeted workers.
template <class F>
inline void parallel_for(ThreadPool *pool, int n, F body) {
  const int W = std::max(1, (int)pool->size());
  const int width = (n + W - 1) / W;
  vector<future<void>> r;
  for (int a = 0; a < n; a += width)
    r.push_back(pool->enqueue([a, b = std::min(a + width, n), &body]() { body(a, b); }));
  joinNclean(r);
}

}  // namespace emp::ag
#endif // EMP_AG_BACKEND_RUNTIME_H__
