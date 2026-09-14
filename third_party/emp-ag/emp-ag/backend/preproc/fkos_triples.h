#pragma once

// Historical CutChoose flow: authenticate random a,b,r with the existing
// COT/aShare provider, hash a's COT correlations into products, authenticate
// the products, then correctness sacrifice and leakage-removal bucketing.
// This is the lean hash-reuse variant, NOT literal FKOS Fig.15 initialization.
// See FKOS.md for the proof scope, checks, and parameter derivation.
#include "emp-ag/wrk.h"
#include "emp-ag/backend/preproc/auth_share_pool.h"
#include "emp-ag/backend/crypto/auth_open.h"
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

namespace emp::ag {
enum class TripleBackend : unsigned { HalfGate = 0, FKOS = 1 };
#ifndef EMP_AG_TRIPLE_BACKEND
#define EMP_AG_TRIPLE_BACKEND 0
#endif
static_assert(EMP_AG_TRIPLE_BACKEND == 0 || EMP_AG_TRIPLE_BACKEND == 1);
inline constexpr unsigned configured_triple_backend_id = EMP_AG_TRIPLE_BACKEND;
inline constexpr TripleBackend configured_triple_backend = static_cast<TripleBackend>(EMP_AG_TRIPLE_BACKEND);
inline constexpr const char *configured_triple_backend_name = EMP_AG_TRIPLE_BACKEND ? "FKOS" : "HalfGate";

struct FkosParameters {
  size_t outputs = 0, candidates = 0, cut = 0;
  unsigned sacrifice = 0, combine = 0;
  long double correctness_error = 0, leakage_error = 0;
};
inline long double fkos_log_choose(size_t n, unsigned k) {
  long double out = 0;
  for (unsigned j = 0; j < k; ++j)
    out += std::log2(static_cast<long double>(n - j)) - std::log2(static_cast<long double>(j + 1));
  return out;
}
// Keep the finite-population factors rather than relaxing to ell^{1-B}.
// Correctness: t/C(T*t,T), t=ell*B, after the random subset check (c>=T).
// Leakage: max_g 2^-g * ell*C(g,B)/C(B*ell,B); max at g=2B-1 or 2B.
// These are bucket-analysis bounds, conditional on the leaky-triple model;
// they are not a UC proof of the hash-reuse candidate generator.
inline FkosParameters fkos_parameters(size_t ell, int ssp) {
  wrk::require(ell >= 2 && ell <= INT_MAX / 4 && ssp >= 1 && ssp <= 120,
               "FKOS parameter range");
  FkosParameters best{};
  const long double budget = std::exp2(-static_cast<long double>(ssp));
  for (unsigned b = 2; b <= unsigned(ssp + 2); ++b) {
    const long double leak = std::exp2(std::log2(static_cast<long double>(ell)) +
        fkos_log_choose(2 * b, b) - 2 * b - fkos_log_choose(ell * b, b));
    if (leak >= budget) continue;
    const size_t t = ell * b;
    for (unsigned tsize = 2; tsize <= unsigned(ssp + 2); ++tsize) {
      const long double correctness = std::exp2(std::log2(static_cast<long double>(t)) -
          fkos_log_choose(t * tsize, tsize));
      if (correctness + leak > budget) continue;
      const size_t cut = std::max<size_t>(tsize, std::ceil(3 * std::log2(static_cast<long double>(t))));
      const size_t n = t * tsize + cut;
      if (!best.candidates || n < best.candidates)
        best = {ell, n, cut, tsize, b, correctness, leak};
      break;
    }
  }
  wrk::require(best.candidates && best.candidates <= INT_MAX / 3, "FKOS candidate count overflow");
  return best;
}
inline size_t fkos_uniform(PRG &prg, size_t bound) {
  wrk::require(bound != 0, "FKOS zero permutation bound");
  const uint64_t n = bound, threshold = -n % n;
  uint64_t x;
  do { prg.random_data_unaligned(&x, sizeof(x)); } while (x < threshold);
  return x % n;
}
inline std::vector<size_t> fkos_permutation(size_t count, block seed) {
  std::vector<size_t> out(count); std::iota(out.begin(), out.end(), 0);
  PRG prg(&seed);
  for (size_t i = count; i > 1; --i) std::swap(out[i - 1], out[fkos_uniform(prg, i)]);
  return out;
}
template<int N> struct BooleanTriple { wrk::AuthShare<N> a, b, c; };

template<int N> class FkosTriples {
 public:
  using Share = wrk::AuthShare<N>;
  using Triple = BooleanTriple<N>;
  static constexpr size_t generation_chunk_size = 16384;
  static constexpr size_t check_padding = 256;
  uint64_t batches = 0, candidates_generated = 0, output_triples = 0;
  uint64_t authentication_cots_per_direction = 0, candidate_chunks = 0;
  FkosParameters last_parameters{};
  struct Metric { std::string name; uint64_t sent = 0, received = 0; double ms = 0; };
  std::vector<Metric> metrics;

  FkosTriples(NetIOMP<N> *io, ThreadPool *pool, int party, block session,
              AuthSharePool<N> &authentication, int ssp)
      : io_(io), pool_(pool), party_(party), session_(session), auth_(authentication), ssp_(ssp) {
    (void)fkos_parameters(128, ssp_);
  }
  std::vector<Triple> generate(size_t count) {
    active();
    try { return generate_impl(count); }
    catch (...) { failed_ = true; throw; }
  }
  void finish() { active(); auth_.maybe_flush_cot_check(); finished_ = true; }
#ifdef EMP_AG_LOCAL_TEST_HOOKS
  std::function<void(const char *, uint8_t *, size_t)> local_test_hook;
#endif

 private:
  NetIOMP<N> *io_;
  ThreadPool *pool_;
  int party_;
  block session_;
  AuthSharePool<N> &auth_;
  int ssp_;
  PRG prg_;
  uint64_t epoch_ = 0;
  bool failed_ = false, finished_ = false;

  void active() const { wrk::require(!failed_ && !finished_, "FKOS session finished or aborted"); }
  void hook(const char *stage, std::vector<uint8_t> &bits) {
#ifdef EMP_AG_LOCAL_TEST_HOOKS
    if (local_test_hook) local_test_hook(stage, bits.data(), bits.size());
#else
    (void)stage; (void)bits;
#endif
  }
  std::pair<uint64_t, uint64_t> traffic() const {
    uint64_t sent = 0, received = 0;
    for (int p = 1; p <= N; ++p) if (p != party_)
      for (auto *ch : {io_->ios[p], io_->ios2[p]}) {
        sent += ch->send_counter; received += ch->recv_counter;
      }
    return {sent, received};
  }
  template<class F> auto measured(const char *name, F f) {
    const auto before = traffic(); const auto start = std::chrono::steady_clock::now();
    auto record = [&] {
      auto it = std::find_if(metrics.begin(), metrics.end(), [&](const Metric &m) { return m.name == name; });
      if (it == metrics.end()) { metrics.push_back({name}); it = metrics.end() - 1; }
      const auto after = traffic();
      it->sent += after.first - before.first; it->received += after.second - before.second;
      it->ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    };
    if constexpr (std::is_void_v<std::invoke_result_t<F>>) { f(); record(); }
    else { auto result = f(); record(); return result; }
  }
  block coin(const char *purpose) {
    return measured("permutation-coins", [&] {
      const block sid = RO("emp-ag:CutChoose:v1", session_).absorb(purpose).absorb(batches).squeeze_block();
      return sampleRandom<N>(io_, &prg_, pool_, party_, sid);
    });
  }
  template<class Send, class Receive> void exchange(Send send, Receive receive) {
    std::vector<std::future<void>> jobs;
    for (int peer : peer_order(party_, N)) {
      jobs.push_back(pool_->enqueue([&, peer] { send(peer); io_->flush(peer); }));
      jobs.push_back(pool_->enqueue([&, peer] { receive(peer); }));
    }
    joinNclean(jobs);
  }
  static Share export_share(const AShareBundle<N> &s) {
    Share out{}; out.bit = own_bit(s);
    for (int j = 0; j < N - 1; ++j) { out.mac(j) = s.mac(j); out.key(j) = s.key(j); }
    return out;
  }
  void add_public(Share &s, bool bit) const {
    if (!bit) return;
    if (party_ == 1) s.bit ^= 1;
    else s.key(peer_slot(party_, 1)) ^= auth_.Delta;
  }
  // Historical lean reuse: hash the very COT blocks that authenticate a.
  // SHA-256/emp-tool Hash with an explicit per-pair, per-candidate domain.
  static uint8_t hash_bit(block scope, size_t index, block value) {
    std::array<unsigned char, 40> input{}; std::memcpy(input.data(), &scope, 16);
    for (unsigned j = 0; j < 8; ++j) input[16 + j] = uint64_t(index) >> (8 * j);
    std::memcpy(input.data() + 24, &value, 16);
    unsigned char digest[Hash::DIGEST_SIZE];
    Hash::hash_once(digest, input.data(), input.size()); return digest[0] & 1;
  }
  std::vector<uint8_t> multiply(const std::vector<Triple> &triples, size_t count, size_t offset) {
    std::array<std::vector<uint8_t>, N + 1> corrections, received, v0, w;
    std::vector<std::future<void>> jobs;
    for (int p : peer_order(party_, N)) {
      corrections[p].resize(count); received[p].resize(count); v0[p].resize(count); w[p].resize(count);
      jobs.push_back(pool_->enqueue([&, p] {
        const int slot = peer_slot(party_, p);
        const block ds = RO("emp-ag:CutChoose:H:v1", session_).absorb(batches)
            .absorb(uint64_t(party_)).absorb(uint64_t(p)).squeeze_block();
        const block dr = RO("emp-ag:CutChoose:H:v1", session_).absorb(batches)
            .absorb(uint64_t(p)).absorb(uint64_t(party_)).squeeze_block();
        for (size_t i = 0; i < count; ++i) {
          const auto &t = triples[offset + i];
          const block k = t.a.key(slot);
          v0[p][i] = hash_bit(ds, offset + i, k);
          corrections[p][i] = v0[p][i] ^ t.b.bit ^ hash_bit(ds, offset + i, k ^ auth_.Delta);
          w[p][i] = hash_bit(dr, offset + i, t.a.mac(slot));
        }
      }));
    }
    joinNclean(jobs);
    exchange([&](int p) { io_->send_bool(p, reinterpret_cast<const bool *>(corrections[p].data()), count); },
             [&](int p) { io_->recv_bool(p, reinterpret_cast<bool *>(received[p].data()), count); });
    std::vector<uint8_t> z(count);
    for (size_t i = 0; i < count; ++i) {
      const uint8_t a = triples[offset + i].a.bit;
      z[i] = a & triples[offset + i].b.bit;
      for (int p = 1; p <= N; ++p) if (p != party_) z[i] ^= v0[p][i] ^ w[p][i] ^ (a & received[p][i]);
    }
    if (offset == 0) hook("products", z);
    hook("all-products", z);
    return z;
  }
  // Raw z is untrusted. Echo-open z+r, authenticate c=r+(z+r), and only
  // then certify correctness by subset checks and bucket sacrifice.
  std::vector<uint8_t> raw_open(const std::vector<uint8_t> &mine) {
    const size_t bytes = (mine.size() + 7) / 8;
    std::vector<uint8_t> packed(bytes);
    for (size_t i = 0; i < mine.size(); ++i) packed[i / 8] |= mine[i] << (i % 8);
    std::array<std::vector<uint8_t>, N + 1> views;
    uint8_t *ptrs[N + 1]{};
    for (int p = 1; p <= N; ++p) { views[p].resize(bytes); ptrs[p] = views[p].data(); }
    EchoBC<N> echo(io_, pool_, party_);
    echo.all_bcast(packed.data(), bytes, ptrs); echo.finalize();
    std::vector<uint8_t> out(mine.size());
    for (int p = 1; p <= N; ++p)
      for (size_t i = 0; i < mine.size(); ++i) out[i] ^= (views[p][i / 8] >> (i % 8)) & 1;
    return out;
  }
  template<class Access> std::vector<uint8_t> open(size_t count, const char *stage, Access access) {
    std::vector<uint8_t> bits(count);
    for (size_t i = 0; i < count; ++i) bits[i] = access(i).bit;
    hook(stage, bits);
    return measured(stage, [&] { return checked_open_bits(io_, pool_, party_, auth_.Delta, session_,
        std::string("FKOS-") + stage, epoch_++, bits,
        [&](int p, size_t i) { return access(i).mac(peer_slot(party_, p)); },
        [&](int p, size_t i) { return access(i).key(peer_slot(party_, p)); }); });
  }
  std::vector<Triple> generate_impl(size_t count) {
    last_parameters = fkos_parameters(count, ssp_); const auto cfg = last_parameters;
    ++batches; candidates_generated += cfg.candidates;
    std::vector<Triple> candidates(cfg.candidates);
    // Small COT-generation chunks, ONE large globally shuffled bucket batch.
    // The 256 discarded aShares mask check2's 128-bit linear disclosure;
    // the random 128x256 padding matrix has rank failure <2^-128. They
    // are never candidates, products, or output bits. All native checks run.
    for (size_t offset = 0; offset < cfg.candidates; offset += generation_chunk_size) {
      const size_t n = std::min(generation_chunk_size, cfg.candidates - offset);
      AShareBundleVec<N> bits;
      measured("cot-authenticate-check", [&] {
        auth_.draw(3 * n + check_padding, bits);
      });
      ++candidate_chunks; authentication_cots_per_direction += 3 * n + check_padding;
      for (size_t i = 0; i < n; ++i)
        candidates[offset + i] = {export_share(bits[i]), export_share(bits[n + i]), export_share(bits[2 * n + i])};
    }
    // Keep Ferret's long-lived stream open across memory chunks. Finish its
    // checks ONCE before any multiplication/hash or product-mask disclosure.
    measured("cot-authenticate-check", [&] { auth_.maybe_flush_cot_check(); });
    for (size_t offset = 0; offset < cfg.candidates; offset += generation_chunk_size) {
      const size_t n = std::min(generation_chunk_size, cfg.candidates - offset);
      auto z = measured("ot-multiply", [&] { return multiply(candidates, n, offset); });
      for (size_t i = 0; i < n; ++i) z[i] ^= candidates[offset + i].c.bit;
      const auto d = measured("authenticate-products", [&] { return raw_open(z); });
      for (size_t i = 0; i < n; ++i) add_public(candidates[offset + i].c, d[i]);
    }
    auto order = fkos_permutation(cfg.candidates, coin("cut"));
    const auto cut = open(3 * cfg.cut, "cut-open", [&](size_t j) -> const Share & {
      const auto &t = candidates[order[j / 3]]; return j % 3 == 0 ? t.a : j % 3 == 1 ? t.b : t.c;
    });
    for (size_t i = 0; i < cfg.cut; ++i)
      wrk::require((cut[3 * i] & cut[3 * i + 1]) == cut[3 * i + 2], "FKOS cut: incorrect authenticated triple");
    order.erase(order.begin(), order.begin() + cfg.cut);

    const auto shuffle = fkos_permutation(order.size(), coin("sacrifice"));
    std::vector<size_t> permutation(order.size());
    for (size_t i = 0; i < order.size(); ++i) permutation[i] = order[shuffle[i]];
    const size_t heads = count * cfg.combine, pairs = heads * (cfg.sacrifice - 1);
    auto pair = [&](size_t j) {
      const size_t h = j / (cfg.sacrifice - 1), k = j % (cfg.sacrifice - 1) + 1;
      return std::pair<size_t, size_t>(permutation[h * cfg.sacrifice], permutation[h * cfg.sacrifice + k]);
    };
    const auto differences = open(2 * pairs, "sacrifice-open", [&](size_t j) {
      const auto [h, k] = pair(j / 2);
      return j % 2 == 0 ? wrk::xor_share(candidates[h].a, candidates[k].a) : wrk::xor_share(candidates[h].b, candidates[k].b);
    });
    const auto checks = open(pairs, "sacrifice-check", [&](size_t j) {
      const auto [h, k] = pair(j);
      Share v = wrk::xor_share(candidates[h].c, candidates[k].c);
      if (differences[2 * j]) v = wrk::xor_share(v, candidates[k].b);
      if (differences[2 * j + 1]) v = wrk::xor_share(v, candidates[k].a);
      add_public(v, differences[2 * j] & differences[2 * j + 1]); return v;
    });
    for (auto v : checks) wrk::require(v == 0, "FKOS sacrifice: incorrect authenticated triple");
    // Compact in increasing source order, safely in-place. A fresh permutation
    // follows, so no particular survivor ordering needs to be preserved.
    order.resize(heads);
    for (size_t i = 0; i < heads; ++i) order[i] = permutation[i * cfg.sacrifice];
    std::sort(order.begin(), order.end());
    for (size_t i = 0; i < heads; ++i) candidates[i] = candidates[order[i]];
    candidates.resize(heads);

    order = fkos_permutation(heads, coin("combine"));
    const size_t joins = count * (cfg.combine - 1);
    const auto dy = open(joins, "combine-open", [&](size_t j) {
      const size_t bucket = j / (cfg.combine - 1), k = j % (cfg.combine - 1) + 1;
      return wrk::xor_share(candidates[order[bucket * cfg.combine]].b, candidates[order[bucket * cfg.combine + k]].b);
    });
    std::vector<Triple> out(count);
    for (size_t i = 0; i < count; ++i) {
      out[i] = candidates[order[i * cfg.combine]];
      for (size_t k = 1; k < cfg.combine; ++k) {
        const auto &t = candidates[order[i * cfg.combine + k]];
        out[i].a = wrk::xor_share(out[i].a, t.a); out[i].c = wrk::xor_share(out[i].c, t.c);
        if (dy[i * (cfg.combine - 1) + k - 1]) out[i].c = wrk::xor_share(out[i].c, t.a);
      }
    }
    output_triples += count; return out;
  }
};
}  // namespace emp::ag
