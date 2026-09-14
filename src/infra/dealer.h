// src/infra/dealer.h — the DEMO trusted dealer: every stand-in for a real offline
// phase lives here, and only here.
//
// A real deployment replaces this file (and nothing else): MASCOT/Overdrive for
// the arithmetic authenticated shares, real malicious COT + an edaBit protocol
// for the Boolean/edaBit correlations, and a distributed KeyGen. The protocol
// layer touches the dealer only through deal_* and my_alpha_* / my_delta.
//
// The fake works by determinism, not secrecy: every party runs the IDENTICAL
// seeded stream, COMPUTES the whole split, and only ever USES its own slot —
// the visibility a real share would give. Correspondingly every party's
// authentication key (arithmetic alpha, Boolean Delta) is dealer-known.
#ifndef MLDSA_DEALER_H
#define MLDSA_DEALER_H

#include "bdoz.h"
#include <emp-ag/wrk.h> // emp::wrk::AuthShare — the GMW/WRK Boolean share type

#include <array>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

namespace mldsa {

// Authenticated arithmetic and Boolean representations of the same value
// (an edaBit). BitWidth/IsUnsigned are part of the type so consumers cannot
// silently mix encodings.
template <int nP, int BitWidth, bool IsUnsigned> struct SharePair {
  static_assert(BitWidth > 0, "SharePair: BitWidth must be positive");
  static constexpr int bit_width = BitWidth;
  static constexpr bool is_unsigned = IsUnsigned;
  std::vector<FqShare<nP>> fq_share;     // pairwise-authenticated F_q share
  std::vector<RingShare<nP>> ring_share; // same value over Z_{2^RING_W}
  emp::wrk::AuthShareVec<nP> two_share;  // GMW/WRK authenticated Boolean bits, LSB first
};

// Ring edaBit (paper: eDaBit_{2^k,k}): a widened R^ in Z_{2^RING_W} and the low
// KB bits of R^ as authenticated Boolean shares.
template <int nP> struct RingEdabits {
  std::vector<RingShare<nP>> ring_share;
  emp::wrk::AuthShareVec<nP> two_share; // KB bits per value, little-endian
};

// The ONLY holder of every party's slopes. In a real BDOZ offline each party
// samples its own alpha_j and never reveals it; the pairwise (key, MAC) pairs
// on a peer's share are produced by correlated OT or homomorphic encryption.
// Replacing this struct is the de-cheat: protocol code only ever touches
// dealer.my_alpha_f / my_alpha_r and the deal_* entry points, all of which a
// real offline can serve without any party holding another's slope.
//
// The fake works by determinism, not secrecy: every party runs the IDENTICAL
// seeded rng, COMPUTES the whole split, and only ever USES its own slot
// [party] — the visibility a real ⟨·⟩ share would give. NOT how real secret
// sharing works; a self-consistent, checkable stand-in.
// Deterministic per-party GMW/WRK authentication key with the aShare pinned-bit
// profile (bit0 = 1, bit1 = party==1 ? nP%2 : 1). Derived from a FIXED master
// block so a FakeDealer can reconstruct EVERY party's Delta and fabricate
// consistent Boolean authenticated shares -- the same demo cheat as the dealer
// knowing every party's arithmetic slope. A real edaBit protocol keeps Delta
// private and produces the share pairs without a dealer.
template <int nP> inline emp::block deterministic_delta(int party) {
  static const emp::block master = emp::makeBlock(0x6d6c647361447fULL, 0x44656c7461ULL);
  emp::PRG prg(&master, party);
  emp::block d;
  prg.random_block(&d, 1);
  std::array<uint8_t, 16> raw{};
  std::memcpy(raw.data(), &d, 16);
  const uint8_t bit1 = (party == 1) ? (uint8_t)(nP % 2) : (uint8_t)1;
  raw[0] = (uint8_t)((raw[0] & ~3U) | 1U | (bit1 << 1));
  std::memcpy(&d, raw.data(), 16);
  return d;
}

template <int nP> struct FakeDealer {
  int party;
  std::mt19937 rng;  // shared seed: every party runs the identical stream
  std::array<RingVal, nP + 1> alpha_r{};  // CHEAT: every party's ring BDOZ slope
  RingVal my_alpha_r;                     // this party's own ring slope
  std::array<FqMac, nP + 1> alpha_f{};    // CHEAT: every party's R F_q slopes
  FqMac my_alpha_f{};                     // this party's own F_q slopes
  std::array<emp::block, nP + 1> deltas{}; // CHEAT: every party's GMW/WRK Delta
  emp::block my_delta{};                   // this party's own Delta

  FakeDealer(int party, uint32_t seed) : party(party), rng(seed) {
    for (int p = 1; p <= nP; ++p)
      alpha_r[p] = ring_rand();
    my_alpha_r = alpha_r[party];
    for (int p = 1; p <= nP; ++p)
      for (int i = 0; i < FQ_MAC_REP; ++i)
        alpha_f[p][i] = rng() % Q;
    my_alpha_f = alpha_f[party];
    for (int p = 1; p <= nP; ++p)
      deltas[p] = deterministic_delta<nP>(p);
    my_delta = deltas[party];
  }

  // A uniform value in [0, range) from the shared stream (identical at every
  // party). Modulo bias is negligible for the demo ranges used here.
  uint64_t draw_uniform(uint64_t range) {
    const uint64_t r = ((uint64_t)rng() << 32) | (uint64_t)rng();
    return range ? r % range : 0;
  }

  // A 128-bit block from the shared stream (identical at every party).
  emp::block block_rand() {
    uint64_t lo = ((uint64_t)rng() << 32) | (uint64_t)rng();
    uint64_t hi = ((uint64_t)rng() << 32) | (uint64_t)rng();
    return emp::makeBlock(hi, lo);
  }

  // Deal an authenticated GMW/WRK Boolean sharing of the low `width` bits of
  // `value` (little-endian). Every ordered pair (holder i, verifier j) gets a
  // one-time key K at P_j and MAC K ^ x_i*Delta_j at P_i, so the receiver of an
  // opening verifies locally -- the same pairwise structure as deal_ring, over
  // GF(2) with Delta instead of the ring slope. Returns THIS party's shares.
  emp::wrk::AuthShareVec<nP> deal_bits(uint64_t value, int width) {
    emp::wrk::AuthShareVec<nP> out((size_t)width);
    for (int k = 0; k < width; ++k) {
      const uint8_t x = (uint8_t)((value >> k) & 1);
      std::array<uint8_t, nP + 1> sh{};
      uint8_t acc = 0;
      for (int p = 1; p < nP; ++p) {
        sh[p] = (uint8_t)(rng() & 1);
        acc ^= sh[p];
      }
      sh[nP] = (uint8_t)(x ^ acc);
      emp::wrk::AuthShare<nP>& a = out[(size_t)k];
      a.bit = sh[party];
      for (int i = 1; i <= nP; ++i)
        for (int j = 1; j <= nP; ++j) {
          if (i == j)
            continue;
          const emp::block K = block_rand();
          if (party == j)
            a.key(emp::wrk::peer_slot(j, i)) = K;
          if (party == i)
            a.mac(emp::wrk::peer_slot(i, j)) = K ^ (emp::select_mask[sh[i]] & deltas[j]);
        }
    }
    return out;
  }

  // BDOZ over the widened ring: additive shares of v, and for every ordered
  // pair (j, i) a key K_{j,i} at P_j and MAC K_{j,i} + alpha_j * v_i at P_i.
  // Every party walks the identical stream and keeps only its own slots.
  RingShare<nP> deal_ring(const RingVal& v) {
    const auto vs = ring_split(v);
    RingShare<nP> out;
    out.val = vs[party];
    for (int j = 1; j <= nP; ++j)
      for (int i = 1; i <= nP; ++i) {
        if (i == j) continue;
        const RingVal k = ring_rand();
        if (party == j) out.key[i] = k;
        if (party == i) out.mac[j] = ring_add(k, ring_mul(alpha_r[j], vs[i]));
      }
    return out;
  }

  // Pairwise-authenticated F_q share of v (same structure as deal_ring).
  FqShare<nP> deal_fq(uint32_t v) {
    const auto vs = arith_split(v);
    FqShare<nP> out;
    out.val = vs[party];
    for (int j = 1; j <= nP; ++j)
      for (int i = 1; i <= nP; ++i) {
        if (i == j) continue;
        FqMac k{};
        for (int r = 0; r < FQ_MAC_REP; ++r)
          k[r] = rng() % Q;
        if (party == j) out.key[i] = k;
        if (party == i)
          for (int r = 0; r < FQ_MAC_REP; ++r)
            out.mac[j][r] = fq_add(k[r], fq_mul(alpha_f[j][r], vs[i]));
      }
    return out;
  }

  // A small-norm secret polynomial for KeyGen: sample each coefficient uniformly
  // in [-eta, eta] and deal its F_q and ring shares. Since the dealer knows the
  // plaintext (it sampled it), it also returns `plain` so KeyGen can compute the
  // public t = A*s + e directly, with no share matrix-multiply and no opening.
  struct SecretPoly {
    std::vector<int32_t> plain;         // centered coefficients in [-eta, eta]
    std::vector<FqShare<nP>> fq;
    std::vector<RingShare<nP>> ring;
  };
  SecretPoly deal_secret_poly(int count, int eta) {
    SecretPoly out;
    out.plain.resize((size_t)count);
    out.fq.resize((size_t)count);
    out.ring.resize((size_t)count);
    for (int i = 0; i < count; ++i) {
      const int32_t v = (int32_t)draw_uniform(2 * (uint64_t)eta + 1) - eta;
      out.plain[(size_t)i] = v;
      out.fq[(size_t)i] = deal_fq((uint32_t)(((v % Q) + Q) % Q));
      out.ring[(size_t)i] = deal_ring(ring_from_i64(v));
    }
    return out;
  }

  // Uniform element of Z_{2^RING_W} from the shared stream (same at every party).
  RingVal ring_rand() {
    RingVal r;
    for (int i = 0; i < RING_WORDS; ++i)
      r.w[i] = ((uint64_t)rng() << 32) | rng();
    ring_reduce(r);
    return r;
  }


  // ---- edaBit families (paper step 1) -----------------------------------------
  // Each value is sampled ALREADY in its required range (no rejection sampling)
  // and dealt in every representation the protocol consumes, at once. Zero
  // communication and zero GMW; the real circuits (C_ad, C_prod, C_post) are the
  // only MPC that runs.

  // y uniform in (-2^(W-1), 2^(W-1)] (ML-DSA's left-open interval). Boolean half
  // is the two's-complement W-bit encoding of v = y - 1; arithmetic halves are y.
  template <int BitWidth> SharePair<nP, BitWidth, false> deal_y_edabit(int count) {
    static_assert(BitWidth > 0 && BitWidth < 32, "deal_y_edabit: BitWidth in [1,31]");
    SharePair<nP, BitWidth, false> out;
    out.fq_share.resize((size_t)count);
    out.ring_share.resize((size_t)count);
    out.two_share.resize((size_t)BitWidth * count);
    for (int i = 0; i < count; ++i) {
      const uint64_t vu = draw_uniform(uint64_t{1} << BitWidth); // v = s (two's complement)
      const int64_t s = (vu >> (BitWidth - 1)) ? (int64_t)vu - (int64_t{1} << BitWidth) : (int64_t)vu;
      const int64_t y = s + 1;
      out.fq_share[(size_t)i] = deal_fq((uint32_t)(((y % Q) + Q) % Q));
      out.ring_share[(size_t)i] = deal_ring(ring_from_i64(y));
      auto bits = deal_bits(vu, BitWidth);
      for (int k = 0; k < BitWidth; ++k)
        out.two_share[(size_t)BitWidth * i + k] = bits[(size_t)k];
    }
    return out;
  }

  // R uniform in [0, Q) with its canonical W-bit Boolean encoding (the A2B mask;
  // ring_share stays empty, it is never consumed).
  template <int BitWidth> SharePair<nP, BitWidth, true> deal_fq_edabit(int count) {
    static_assert(BitWidth > 0 && BitWidth < 32, "deal_fq_edabit: BitWidth in [1,31]");
    SharePair<nP, BitWidth, true> out;
    out.fq_share.resize((size_t)count);
    out.two_share.resize((size_t)BitWidth * count);
    for (int i = 0; i < count; ++i) {
      const uint32_t r = (uint32_t)draw_uniform((uint64_t)Q);
      out.fq_share[(size_t)i] = deal_fq(r);
      auto bits = deal_bits(r, BitWidth);
      for (int k = 0; k < BitWidth; ++k)
        out.two_share[(size_t)BitWidth * i + k] = bits[(size_t)k];
    }
    return out;
  }

  // Widened R^ uniform in Z_{2^RING_W}; the low KB bits are the Boolean half.
  template <int KB> RingEdabits<nP> deal_ring_edabit(int count) {
    static_assert(KB > 0 && KB < 32, "deal_ring_edabit: KB in [1,31]");
    constexpr uint64_t KMASK = (uint64_t{1} << KB) - 1;
    RingEdabits<nP> out;
    out.ring_share.resize((size_t)count);
    out.two_share.resize((size_t)KB * count);
    for (int i = 0; i < count; ++i) {
      const RingVal v = ring_rand();
      out.ring_share[(size_t)i] = deal_ring(v);
      auto bits = deal_bits(v.w[0] & KMASK, KB);
      for (int k = 0; k < KB; ++k)
        out.two_share[(size_t)KB * i + k] = bits[(size_t)k];
    }
    return out;
  }

private:
  std::array<RingVal, nP + 1> ring_split(const RingVal& v) {
    std::array<RingVal, nP + 1> s{};
    RingVal acc;
    for (int p = 1; p < nP; ++p) {
      s[p] = ring_rand();
      acc = ring_add(acc, s[p]);
    }
    s[nP] = ring_sub(v, acc);
    return s;
  }

  // Additive split of v over F_Q; slot p is party p's share.
  std::array<uint32_t, nP + 1> arith_split(uint32_t v) {
    std::array<uint32_t, nP + 1> s{};
    uint32_t acc = 0;
    for (int p = 1; p < nP; ++p) {
      s[p] = rng() % Q;
      acc = fq_add(acc, s[p]);
    }
    s[nP] = fq_sub(v, acc); // last share absorbs the rest
    return s;
  }
};

} // namespace mldsa
#endif // MLDSA_DEALER_H
