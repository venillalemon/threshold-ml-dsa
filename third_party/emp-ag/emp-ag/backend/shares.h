#ifndef EMP_AG_BACKEND_SHARES_H__
#define EMP_AG_BACKEND_SHARES_H__
//
// The authenticated-share data model: the physical per-wire share (SecureShare and its
// ADL StorageTraits accessors), the agmpc aliases (AShareBundle + peer_slot), and the
// per-circuit wire container (SecureWires). The section comments below carry the detail.
//
#include "emp-tool/runtime/core/block.h"
#include "emp-tool/runtime/core/block_vector.h"
#include <vector>

namespace emp::ag {

//
// Per-wire authenticated-share storage — the single physical share type the whole stack
// carries.
//
// Layout (slot k indexes peers in increasing party-id order, skipping `self`):
//   km[2*k]     = mac for the k-th peer = M_k[x] = K_k[x] ^ x*Delta_k  (bit0 = x)
//   km[2*k + 1] = key for the k-th peer = K_me[x_k]
// Size = 2*(nP-1) blocks. At nP=2 (K=1) this is exactly 32B = {mac, key}, and the slot-0
// accessors + the XOR-over-one-peer traits below fold to scalar loads (the per-peer loop
// has one iteration; `zero_block ^ x` folds to `x`). That is what keeps the two-party hot
// loop's scalar codegen with no n-party machinery forced onto it (the StorageTraits seam:
// the ADL accessors below -- peer_key_sum / peer_mac_sum / own_bit / peer_mac / own_dot).
//
template <int nP> struct SecureShare {
  static_assert(nP >= 2, "at least 2 parties");
  static constexpr int K = nP - 1;          // peers (everyone but self)
  static constexpr int N = 2 * (nP - 1);    // blocks: (mac,key) per peer
  block km[N];

  block &mac(int slot) { return km[2 * slot]; }
  block &key(int slot) { return km[2 * slot + 1]; }
  const block &mac(int slot) const { return km[2 * slot]; }
  const block &key(int slot) const { return km[2 * slot + 1]; }
};

// StorageTraits accessors, found by ADL from the garble / output-check kernels. At nP=2
// (K=1) the loop is one iteration and `zero_block ^ key(0)` folds to `key(0)`.
template <int nP> inline block peer_key_sum(const SecureShare<nP> &b) {
  block s = zero_block;
  for (int k = 0; k < SecureShare<nP>::K; ++k) s = s ^ b.key(k);
  return s;
}
// XOR of the local MACs over all peers = XOR_k M_k[x^me]. The MAC-side counterpart of
// peer_key_sum, used by the BDOZ->SPDZ collapse (collapse_cgamma_witness, crypto/spdz.h): the SPDZ MAC share of
// x under the shared global key is own_dot(x,Delta_me) ^ peer_mac_sum(x) ^ peer_key_sum(x).
// At nP=2 (K=1) this folds to the single scalar mac(0).
template <int nP> inline block peer_mac_sum(const SecureShare<nP> &b) {
  block s = zero_block;
  for (int k = 0; k < SecureShare<nP>::K; ++k) s = s ^ b.mac(k);
  return s;
}
template <int nP> inline bool own_bit(const SecureShare<nP> &b) {
  return getLSB(b.mac(0));
}
// own_bit * Delta, branchless via select_mask (folds the bit to 0 / Delta), as the
// half-gate and c_gamma kernels use it. Resolves own_bit by ADL on B.
template <class B>
inline block own_dot(const B &b, const block &Delta) {
  return select_mask[own_bit(b)] & Delta;
}
// MAC under a specific peer's key (slot s) -- the shared t_gamma garbler witness
// over the evaluator's slot. At nP=2 the only slot is 0.
template <int nP> inline block peer_mac(const SecureShare<nP> &b, int s) {
  return b.mac(s);
}

template <int nP>
using SecureShareVec =
    std::vector<SecureShare<nP>, default_init_allocator<SecureShare<nP>>>;

// The per-wire share IS emp::ag::SecureShare<nP>. The
// share's ADL traits (peer_key_sum/own_bit/peer_mac) live above
// and are found by ADL from the shared kernels; the peer_slot index helper below is
// agmpc-local (plain ints, no ADL).
template <int nP> using AShareBundle = emp::ag::SecureShare<nP>;

template <int nP> using AShareBundleVec = emp::ag::SecureShareVec<nP>;

// Free-XOR share propagation: recompute a fabric wire's share bundle as the
// componentwise XOR of its (live) inputs. bit0 (the share-bit) rides along.
template <int nP>
inline void xor_share(AShareBundle<nP> &out, const AShareBundle<nP> &a,
                      const AShareBundle<nP> &b) {
  for (int t = 0; t < AShareBundle<nP>::N; ++t) out.km[t] = a.km[t] ^ b.km[t];
}
template <int nP>
inline void zero_share(AShareBundle<nP> &out) {
  for (int t = 0; t < AShareBundle<nP>::N; ++t) out.km[t] = zero_block;
}

// Map peer party-id (1-based) to a slot index in [0, nP-1) within a bundle
// owned by `self`. Caller must ensure peer != self. (agmpc-only; ints, no ADL.)
inline int peer_slot(int self, int peer) {
  return (peer < self) ? (peer - 1) : (peer - 2);
}

// Bit 1 (second-lowest) of a block. Reserved for the half-gate Λ_γ recovery
// invariant LSB1(⊕_p Δ_p) = 1; bit 0 is reserved for share-value encoding
// (SecureShare, above).
inline uint8_t LSB1(const block &b) { return getBit(b, 1); }

// Per-wire authenticated-share + half-gate state. Notation (Pi = local):
//   wire_bundle[w]   = (M_j[λ_w^i], K_i[λ_w^j]) for j ≠ i, AoS-by-wire.
//                      Slot k = peer_slot(party, j) maps peer-id j → slot.
//                      mac(k) / key(k) accessors hide the layout. At nP=3
//                      one bundle = 64B = 1 cache line, so per-wire access
//                      from a topological garble loop hits a single line.
//   λ_w^i            = getLSB(wire_bundle[w].mac(0))  — implicit, since the
//                      bit-0 invariant pins bit0(M)=x and bit0(K)=0 across
//                      all peer slots.
//   Lambda[w]        = Λ_w             (publicly opened mask; populated
//                                       at every party after process_input
//                                       or compute)
//   label0[w]        = m_{w, 0}^i      only at Pi (i ≥ 2); empty at P1
//   eval_label[j][w] = m_{w, Λ_w}^j    for j ≥ 2, only at P1; empty at Pi (i ≥ 2)
//
// MAC relation: M_j[λ_w^i] = K_j[λ_w^i] ⊕ λ_w^i · Δ_j; on Pi the LHS is
// in bundle.mac(slot of j) and on Pj the RHS-key is in (Pj's local)
// bundle.key(slot of i).
template <int nP> struct SecureWires {
  std::vector<unsigned char> Lambda;
  AShareBundleVec<nP> wire_bundle;
  BlockVec label0;
  BlockVec eval_label[nP + 1];

  size_t size() const { return Lambda.size(); }

  // Extract wires [lo, hi) into a fresh bundle. Per-party label vectors
  // (eval_label) are sliced when populated and left empty otherwise,
  // mirroring the layout invariant from process_input / compute.
  SecureWires<nP> slice(size_t lo, size_t hi) const {
    SecureWires<nP> r;
    r.Lambda.assign(Lambda.begin() + lo, Lambda.begin() + hi);
    if (!wire_bundle.empty())
      r.wire_bundle.assign(wire_bundle.begin() + lo, wire_bundle.begin() + hi);
    for (int j = 1; j <= nP; ++j) {
      if (!eval_label[j].empty())
        r.eval_label[j].assign(eval_label[j].begin() + lo,
                               eval_label[j].begin() + hi);
    }
    if (!label0.empty())
      r.label0.assign(label0.begin() + lo, label0.begin() + hi);
    return r;
  }

  // In-place concatenation: append b's wires onto `*this` without allocating a copy.
  void append(const SecureWires<nP> &b) {
    Lambda.insert(Lambda.end(), b.Lambda.begin(), b.Lambda.end());
    wire_bundle.insert(wire_bundle.end(), b.wire_bundle.begin(),
                       b.wire_bundle.end());
    for (int j = 1; j <= nP; ++j) {
      eval_label[j].insert(eval_label[j].end(), b.eval_label[j].begin(),
                           b.eval_label[j].end());
    }
    label0.insert(label0.end(), b.label0.begin(), b.label0.end());
  }
};

}  // namespace emp::ag
#endif  // EMP_AG_BACKEND_SHARES_H__
