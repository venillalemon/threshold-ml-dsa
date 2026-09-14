#ifndef EMP_AG_BACKEND_PROFILING_H__
#define EMP_AG_BACKEND_PROFILING_H__

// Opt-in profiling for the AG-MPC protocol stack.
//
// Compile with -DAG_PROFILE (or define it before including protocol headers)
// to turn on phase timing, triple-pool subphase timing, communication counters,
// and the shared sub-protocol byte counters. In normal builds the macros compile
// away.

#ifdef AG_PROFILE
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <sys/resource.h>

namespace emp::ag {

// This party's send+recv bytes spent in COT extension and half-gate phi
// exchange. Producers increment these counters; draw_and_seed (passes/draw_beaver.h)
// samples them to split aggregate communication into protocol subcomponents.
inline uint64_t g_ag_cot_bytes = 0;
inline uint64_t g_ag_phi_bytes = 0;

// Sub-phase wall accumulators for the layered leaky-AND triple generation
// (compute_inplace resets them per call and prints the [ag-tp] breakdown at
// the end):
//   cot    = COT / aShare draw (abit.compute -- includes check2's folds at nP>2)
//   hg     = half-gate phi exchange + hash (produce_leaky_ands_halfgate)
//   sopen  = leaky-AND s-bit open
//   acheck = d-combine + the LaAND consistency accumulation (the striped F_eq
//            window fold at nP==2 / the windowed-RLC stash at nP>=3) -- the
//            d-combine dominates this bucket now that both are pool-split
//   bkt    = bucket_one_layer (cyclic-shift fold + its d-open)
//   open   = final Beaver x/y open + sigma assembly
inline uint64_t g_ag_tp_cot_ns = 0, g_ag_tp_hg_ns = 0;
inline uint64_t g_ag_tp_sopen_ns = 0, g_ag_tp_acheck_ns = 0;
inline uint64_t g_ag_tp_bkt_ns = 0, g_ag_tp_open_ns = 0;
inline uint64_t agmpc_now_ns() {
  return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch()).count();
}

inline long agmpc_peak_rss_kib() {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
  return ru.ru_maxrss / 1024;
#else
  return ru.ru_maxrss;
#endif
}

}  // namespace emp::ag

// Requires NetIOMP<nP> *io and int party in the call-site scope.
#define AG_PHASE_BEGIN()                                                     \
  auto _agmpc_t = std::chrono::steady_clock::now();                             \
  int64_t _agmpc_c = io->count()
#define AG_PHASE(name)                                                       \
  do {                                                                          \
    auto _n = std::chrono::steady_clock::now();                                 \
    int64_t _c = io->count();                                                   \
    if (party == 1)                                                             \
      std::printf("[ag] %-22s %9.3f ms  %12lld B  peakRSS %8ld KiB\n",       \
                  (name),                                                       \
                  std::chrono::duration<double, std::milli>(_n - _agmpc_t)      \
                      .count(),                                                 \
                  (long long)(_c - _agmpc_c), emp::ag::agmpc_peak_rss_kib());   \
    _agmpc_t = _n;                                                              \
    _agmpc_c = _c;                                                              \
  } while (0)
#else
#define AG_PHASE_BEGIN() ((void)0)
#define AG_PHASE(name) ((void)0)
#endif

#endif  // EMP_AG_BACKEND_PROFILING_H__
