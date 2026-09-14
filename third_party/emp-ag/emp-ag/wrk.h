#pragma once

// Classical multiparty WRK authenticated garbling, Figures 2--3 of
// Wang--Ranellucci--Katz, Global-Scale Secure Multiparty Computation,
// https://eprint.iacr.org/2017/189 (CCS 2017).
//
// Four rows per AND and garbler. Each decrypted row contains a full-width
// authenticated share of its masked output and shares of all output labels.
// The evaluator checks every selected row locally. This is NOT KRRW half-gate
// evaluation with its interactive c_gamma check removed.
//
// This module assumes correctly authenticated preprocessing (including the
// mask products). prepare_local is party-local and performs no networking.
// make_test_material is a centralized TRUSTED-DEALER TEST FIXTURE ONLY.
#include <emp-tool/emp-tool.h>
#include <emp-tool/ir/program.h>
#include <emp-tool/ir/validate.h>
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cstring>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

namespace emp::wrk {

inline void require(bool good, const char *message) {
  if (!good) throw std::runtime_error(message);
}
inline int peer_slot(int self, int peer) {
  return peer < self ? peer - 1 : peer - 2;
}

// Explicit bit: neither a MAC nor Delta reserves low bits. All 128 MAC/key
// bits are used. On holder i, mac(j)=K_j[x_i] XOR x_i Delta_j; key(j)
// authenticates peer j's share under this party's private Delta_i.
template <int nP> struct AuthShare {
  static_assert(nP >= 2);
  uint8_t bit = 0;
  std::array<block, nP - 1> macs{}, keys{};
  block &mac(int slot) { return macs[slot]; }
  block &key(int slot) { return keys[slot]; }
  const block &mac(int slot) const { return macs[slot]; }
  const block &key(int slot) const { return keys[slot]; }
};
template <int nP> using ShareVec = std::vector<AuthShare<nP>>;
template <int nP> using AuthShareVec = ShareVec<nP>;
template <int nP> inline bool own_bit(const AuthShare<nP> &s) { return s.bit; }
template <int nP> AuthShare<nP> xor_share(const AuthShare<nP> &a,
                                        const AuthShare<nP> &b) {
  AuthShare<nP> out;
  out.bit = a.bit ^ b.bit;
  for (int j = 0; j < nP - 1; ++j) {
    out.mac(j) = a.mac(j) ^ b.mac(j);
    out.key(j) = a.key(j) ^ b.key(j);
  }
  return out;
}

template <int nP> struct GcState {
  uint32_t party = 0, fixed_inputs = 0, num_inputs = 0;
  uint32_t num_wires = 0, num_ands = 0;
  block Delta = zero_block, hash_seed = zero_block;
  std::vector<uint8_t> input_masked, late_public_masks;
  ShareVec<nP> late_input_masks;
  BlockVec input_label0;
  std::array<BlockVec, nP + 1> fixed_input_labels;
  ShareVec<nP> output_masks;
  // P1: complete table; a freshly prepared garbler: its own four rows/gate.
  // After offline delivery a garbler's rows may be erased.
  std::vector<uint8_t> rows;
  // P1's authenticated affine bases for all four rows, three shares/AND.
  ShareVec<nP> ev_base, ev_alpha, ev_beta;
};

template <int nP> constexpr size_t row_bytes() { return 1 + size_t(16) * nP; }
template <int nP> size_t row_offset(size_t gate, int garbler, unsigned row) {
  return ((gate * (nP - 1) + size_t(garbler - 2)) * 4 + row) * row_bytes<nP>();
}
inline uint32_t count_gc_ands(const circuit::BooleanProgram &program) {
  return static_cast<uint32_t>(std::count_if(program.gates.begin(),
      program.gates.end(), [](const auto &g) { return g.op == circuit::Op::And; }));
}

namespace detail {
inline void put_block(uint8_t *dst, block b) { std::memcpy(dst, &b, 16); }
inline block get_block(const uint8_t *src) {
  block b; std::memcpy(&b, src, 16); return b;
}
inline void put_u64(uint8_t *dst, uint64_t x) {
  for (int i = 0; i < 8; ++i) dst[i] = uint8_t(x >> (8 * i));
}

// H(label_alpha, label_beta, gate, row) in the paper, instantiated as a
// domain-separated SHAKE256 XOF. Bind the execution seed, exact output length,
// garbler role, gate index and row, and preserve the two label positions.
// No AES correlation-robustness assumption or speculative row aggregation.
inline void row_pad(std::span<uint8_t> out, block seed, uint64_t gate,
                    uint64_t row, uint64_t garbler, block a, block b) {
  static constexpr char domain[] = "ThresholdMLDSA/WRK17b/four-row/SHAKE256/v1";
  std::array<uint8_t, 80> input{};
  put_block(input.data(), seed);
  put_u64(input.data() + 16, gate);
  put_u64(input.data() + 24, row);
  put_u64(input.data() + 32, garbler);
  put_u64(input.data() + 40, out.size());
  put_block(input.data() + 48, a);
  put_block(input.data() + 64, b);
  using Handle = std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)>;
  Handle ctx(EVP_MD_CTX_new(), EVP_MD_CTX_free);
  require(ctx && EVP_DigestInit_ex(ctx.get(), EVP_shake256(), nullptr) == 1 &&
      EVP_DigestUpdate(ctx.get(), domain, sizeof(domain)) == 1 &&
      EVP_DigestUpdate(ctx.get(), input.data(), input.size()) == 1 &&
      EVP_DigestFinalXOF(ctx.get(), out.data(), out.size()) == 1,
      "WRK SHAKE256 failed");
}

template <int nP>
AuthShare<nP> row_share(const AuthShare<nP> &base,
                        const AuthShare<nP> &alpha,
                        const AuthShare<nP> &beta,
                        bool u, bool v, int party, block Delta) {
  auto s = base;
  if (u) s = xor_share(s, beta);
  if (v) s = xor_share(s, alpha);
  // Public uv is injected into holder 1 only (Fig.2, step4b/c). The holder's
  // MACs stay unchanged; each verifier changes the key for that holder.
  if (u && v) {
    if (party == 1) s.bit ^= 1;
    else s.key(peer_slot(party, 1)) = s.key(peer_slot(party, 1)) ^ Delta;
  }
  return s;
}
}  // namespace detail

// Fully party-local preprocessing interface. input_masks / and_masks / products
// must come from a malicious-secure authenticated preprocessing mechanism under
// this same party's Delta. products[a] authenticates lambda_alpha*lambda_beta,
// NOT a fresh independent Beaver triple's c. No global plaintext is accepted.
template <int nP> struct LocalPreprocessing {
  int party = 0;
  block Delta = zero_block;
  ShareVec<nP> input_masks, and_masks, products;
  BlockVec input_labels0, and_labels0; // sampled by each garbler independently
};

template <int nP>
GcState<nP> prepare_local(const circuit::BooleanProgram &program,
                         const LocalPreprocessing<nP> &material,
                         block hash_seed, uint32_t fixed_inputs,
                         std::span<const uint8_t> fixed_masked,
                         std::span<const uint8_t> late_public_masks) {
  circuit::validate_program(program);
  const int party = material.party;
  const uint32_t na = count_gc_ands(program);
  require(party >= 1 && party <= nP && fixed_inputs <= program.num_inputs &&
      fixed_masked.size() == fixed_inputs &&
      late_public_masks.size() == program.num_inputs - fixed_inputs,
      "WRK preprocessing input shape mismatch");
  require(program.wire_reuse == circuit::WireReuse::None &&
      program.num_wires <= INT_MAX && program.gates.size() <= INT_MAX,
      "WRK requires a dense, non-recycled circuit layout");
  require(material.input_masks.size() == program.num_inputs &&
      material.and_masks.size() == na && material.products.size() == na,
      "WRK authenticated preprocessing shape mismatch");
  require((party == 1 && material.input_labels0.empty() && material.and_labels0.empty()) ||
      (party != 1 && material.input_labels0.size() == program.num_inputs &&
       material.and_labels0.size() == na), "WRK garbler label shape mismatch");
  require(!cmpBlock(&material.Delta, &zero_block, 1), "WRK zero authentication key");
  for (auto bit : fixed_masked) require(bit <= 1, "WRK non-bit fixed mask");
  for (auto bit : late_public_masks) require(bit <= 1, "WRK non-bit late mask");
  for (auto *shares : {&material.input_masks, &material.and_masks, &material.products})
    for (const auto &share : *shares) require(share.bit <= 1, "WRK non-bit authenticated share");
  GcState<nP> s;
  s.party = party; s.fixed_inputs = fixed_inputs; s.num_inputs = program.num_inputs;
  s.num_wires = program.num_wires; s.num_ands = na;
  s.Delta = material.Delta; s.hash_seed = hash_seed;
  s.input_masked.assign(program.num_inputs, 0);
  std::copy(fixed_masked.begin(), fixed_masked.end(), s.input_masked.begin());
  s.late_public_masks.assign(late_public_masks.begin(), late_public_masks.end());
  s.late_input_masks.assign(material.input_masks.begin() + fixed_inputs,
                            material.input_masks.end());
  s.input_label0 = material.input_labels0;
  ShareVec<nP> wires(program.num_wires);
  std::copy(material.input_masks.begin(), material.input_masks.end(), wires.begin());
  BlockVec labels;
  if (party == 1) {
    s.ev_base.resize(na); s.ev_alpha.resize(na); s.ev_beta.resize(na);
  } else {
    labels.resize(program.num_wires);
    std::copy(material.input_labels0.begin(), material.input_labels0.end(), labels.begin());
    s.rows.resize(size_t(4) * na * row_bytes<nP>());
  }
  size_t ai = 0;
  for (const auto &g : program.gates) {
    using circuit::Op;
    if (g.op == Op::Xor) {
      wires[g.out] = xor_share(wires[g.in0], wires[g.in1]);
      if (party != 1) labels[g.out] = labels[g.in0] ^ labels[g.in1];
    } else if (g.op == Op::Not) {
      wires[g.out] = wires[g.in0];
      if (party != 1) labels[g.out] = labels[g.in0] ^ s.Delta;
    } else if (g.op == Op::Const0 || g.op == Op::Const1) {
      wires[g.out] = AuthShare<nP>{};
      if (party != 1) labels[g.out] = g.op == Op::Const1 ? s.Delta : zero_block;
    } else if (g.op == Op::And) {
      // Save operands before writing output; the circuit ABI forbids AND
      // outputs from aliasing live operands.
      const auto alpha = wires[g.in0], beta = wires[g.in1];
      wires[g.out] = material.and_masks[ai];
      const auto base = xor_share(material.products[ai], wires[g.out]);
      if (party == 1) {
        s.ev_base[ai] = base; s.ev_alpha[ai] = alpha; s.ev_beta[ai] = beta;
      } else {
        labels[g.out] = material.and_labels0[ai];
        for (unsigned row = 0; row < 4; ++row) {
          const bool u = row >> 1, v = row & 1;
          const auto share = detail::row_share(base, alpha, beta, u, v, party, s.Delta);
          const size_t at = (4 * ai + row) * row_bytes<nP>();
          auto target = std::span(s.rows).subspan(at, row_bytes<nP>());
          detail::row_pad(target, hash_seed, ai, row, party,
              labels[g.in0] ^ (select_mask[u] & s.Delta),
              labels[g.in1] ^ (select_mask[v] & s.Delta));
          std::array<uint8_t, row_bytes<nP>()> plaintext{};
          plaintext[0] = share.bit;
          block self = labels[g.out] ^ (select_mask[share.bit] & s.Delta);
          for (int peer = 1; peer <= nP; ++peer) if (peer != party) {
            const int slot = peer_slot(party, peer);
            detail::put_block(plaintext.data() + 1 + 16 * slot, share.mac(slot));
            self = self ^ share.key(slot);
          }
          detail::put_block(plaintext.data() + 1 + 16 * (nP - 1), self);
          for (size_t j = 0; j < target.size(); ++j) target[j] ^= plaintext[j];
        }
      }
      ++ai;
    } else throw std::runtime_error("WRK unsupported gate");
  }
  for (auto w : program.outputs) s.output_masks.push_back(wires[w]);
  return s;
}

template <int nP>
void bind_public_inputs(GcState<nP> &s, std::span<const uint8_t> bits) {
  require(s.fixed_inputs <= s.num_inputs && s.input_masked.size() == s.num_inputs &&
      bits.size() == s.num_inputs - s.fixed_inputs &&
      s.late_public_masks.size() == bits.size(), "WRK public input width mismatch");
  for (size_t i = 0; i < bits.size(); ++i) {
    require(bits[i] <= 1 && s.late_public_masks[i] <= 1, "WRK non-bit public input");
    s.input_masked[s.fixed_inputs + i] = bits[i] ^ s.late_public_masks[i];
  }
}
template <int nP>
BlockVec selected_late_labels(GcState<nP> &s, std::span<const uint8_t> bits) {
  require(s.party >= 2 && s.party <= nP && s.input_label0.size() == s.num_inputs,
           "WRK invalid garbler state");
  bind_public_inputs(s, bits);
  BlockVec out(bits.size());
  for (size_t i = 0; i < bits.size(); ++i) {
    const size_t w = s.fixed_inputs + i;
    out[i] = s.input_label0[w] ^ (select_mask[s.input_masked[w]] & s.Delta);
  }
  return out;
}

template <int nP>
std::vector<uint8_t> evaluate(GcState<nP> &s, const circuit::BooleanProgram &program,
                             const std::array<BlockVec, nP + 1> &input_labels,
                             [[maybe_unused]] ThreadPool *pool = nullptr) {
  circuit::validate_program(program);
  require(program.wire_reuse == circuit::WireReuse::None &&
      s.party == 1 && s.num_inputs == program.num_inputs &&
      s.num_wires == program.num_wires && s.num_ands == count_gc_ands(program) &&
      s.input_masked.size() == s.num_inputs &&
      s.ev_base.size() == s.num_ands && s.ev_alpha.size() == s.num_ands &&
      s.ev_beta.size() == s.num_ands &&
      s.rows.size() == size_t(4) * (nP - 1) * s.num_ands * row_bytes<nP>(),
      "WRK evaluator state shape mismatch");
  for (auto *shares : {&s.ev_base, &s.ev_alpha, &s.ev_beta})
    for (const auto &share : *shares) require(share.bit <= 1, "WRK non-bit evaluator share");
  std::vector<uint8_t> masked(s.num_wires);
  std::copy(s.input_masked.begin(), s.input_masked.end(), masked.begin());
  for (auto x : s.input_masked) require(x <= 1, "WRK non-bit masked input");
  std::array<BlockVec, nP + 1> labels;
  for (int p = 2; p <= nP; ++p) {
    require(input_labels[p].size() == s.num_inputs, "WRK input label width mismatch");
    labels[p].resize(s.num_wires);
    std::copy(input_labels[p].begin(), input_labels[p].end(), labels[p].begin());
  }
  size_t ai = 0;
  // One row from each garbler; plaintext rows never leave the evaluator.
  std::array<std::array<uint8_t, row_bytes<nP>()>, nP + 1> plain{};
  for (const auto &g : program.gates) {
    using circuit::Op;
    if (g.op == Op::And) {
      const bool u = masked[g.in0], v = masked[g.in1];
      const unsigned row = 2 * unsigned(u) + unsigned(v);
      const auto mine = detail::row_share(s.ev_base[ai], s.ev_alpha[ai],
                                          s.ev_beta[ai], u, v, 1, s.Delta);
      uint8_t output_masked = mine.bit;
      for (int p = 2; p <= nP; ++p) {
        detail::row_pad(plain[p], s.hash_seed, ai, row, p,
                        labels[p][g.in0], labels[p][g.in1]);
        const size_t at = row_offset<nP>(ai, p, row);
        for (size_t j = 0; j < row_bytes<nP>(); ++j) plain[p][j] ^= s.rows[at + j];
        require(plain[p][0] <= 1, "WRK selected-row noncanonical bit");
        const int slot = peer_slot(1, p);
        const block expected = mine.key(slot) ^ (select_mask[plain[p][0]] & s.Delta);
        const block supplied = detail::get_block(plain[p].data() + 1 + 16 * peer_slot(p, 1));
        require(cmpBlock(&expected, &supplied, 1), "WRK selected-row MAC verification failed");
        output_masked ^= plain[p][0];
      }
      // All selected-row mask shares are authenticated before they affect
      // the public mask fabric. Reconstruct each garbler's selected label.
      masked[g.out] = output_masked;
      for (int dest = 2; dest <= nP; ++dest) {
        block value = detail::get_block(plain[dest].data() + 1 + 16 * (nP - 1));
        value = value ^ mine.mac(peer_slot(1, dest));
        for (int source = 2; source <= nP; ++source) if (source != dest)
          value = value ^ detail::get_block(plain[source].data() + 1 +
                                             16 * peer_slot(source, dest));
        labels[dest][g.out] = value;
      }
      ++ai;
    } else if (g.op == Op::Xor) {
      masked[g.out] = masked[g.in0] ^ masked[g.in1];
      for (int p = 2; p <= nP; ++p) labels[p][g.out] = labels[p][g.in0] ^ labels[p][g.in1];
    } else if (g.op == Op::Not) {
      masked[g.out] = masked[g.in0] ^ 1;
      for (int p = 2; p <= nP; ++p) labels[p][g.out] = labels[p][g.in0];
    } else if (g.op == Op::Const0 || g.op == Op::Const1) {
      masked[g.out] = g.op == Op::Const1;
      for (int p = 2; p <= nP; ++p) labels[p][g.out] = zero_block;
    } else throw std::runtime_error("WRK unsupported gate");
  }
  std::vector<uint8_t> outputs;
  for (auto w : program.outputs) outputs.push_back(masked[w]);
  // These are still masked outputs. Decode with authenticated output masks
  // only after this function succeeds; neither gate checks nor any function
  // in this module require evaluator-to-party communication.
  return outputs;
}

namespace test_dealer_detail {
template <int nP>
std::array<AuthShare<nP>, nP + 1> share_bit(
    uint8_t clear, const std::array<LocalPreprocessing<nP>, nP + 1> &materials,
    PRG &prg) {
  std::array<AuthShare<nP>, nP + 1> result{};
  uint8_t remaining = clear;
  for (int p = 1; p < nP; ++p) {
    prg.random_data_unaligned(&result[p].bit, 1);
    result[p].bit &= 1; remaining ^= result[p].bit;
  }
  result[nP].bit = remaining;
  for (int holder = 1; holder <= nP; ++holder)
    for (int verifier = 1; verifier <= nP; ++verifier) if (holder != verifier) {
      block key; prg.random_block(&key, 1);
      result[verifier].key(peer_slot(verifier, holder)) = key;
      result[holder].mac(peer_slot(holder, verifier)) =
          key ^ (select_mask[result[holder].bit] & materials[verifier].Delta);
    }
  return result;
}
}  // namespace test_dealer_detail

// TEST ONLY: this centralized dealer knows every input, share and private key.
// Its per-party files are useful for online integration testing, not a secure
// distributed preprocessing realization. Production may call prepare_local
// with individually obtained, authenticated correlated material instead.
template <int nP>
std::array<GcState<nP>, nP + 1> make_test_material(
    const circuit::BooleanProgram &program, std::span<const uint8_t> fixed_bits,
    uint32_t fixed_inputs) {
  circuit::validate_program(program);
  require(program.wire_reuse == circuit::WireReuse::None,
           "WRK test dealer requires dense non-recycled wires");
  require(fixed_inputs <= program.num_inputs && fixed_bits.size() == fixed_inputs,
           "WRK fixed input width mismatch");
  std::array<LocalPreprocessing<nP>, nP + 1> materials;
  const uint32_t na = count_gc_ands(program);
  PRG prg; block seed; prg.random_block(&seed, 1);
  for (int p = 1; p <= nP; ++p) {
    auto &m = materials[p]; m.party = p;
    do { prg.random_block(&m.Delta, 1); }
    while (cmpBlock(&m.Delta, &zero_block, 1));
    m.input_masks.resize(program.num_inputs); m.and_masks.resize(na); m.products.resize(na);
    if (p != 1) {
      m.input_labels0.resize(program.num_inputs); m.and_labels0.resize(na);
      prg.random_block(m.input_labels0.data(), m.input_labels0.size());
      prg.random_block(m.and_labels0.data(), m.and_labels0.size());
    }
  }
  std::vector<uint8_t> masks(program.num_wires), fixed_masked(fixed_inputs), late_masks;
  for (uint32_t w = 0; w < program.num_inputs; ++w) {
    prg.random_data_unaligned(&masks[w], 1); masks[w] &= 1;
    const auto share = test_dealer_detail::share_bit<nP>(masks[w], materials, prg);
    for (int p = 1; p <= nP; ++p) materials[p].input_masks[w] = share[p];
    if (w < fixed_inputs) {
      require(fixed_bits[w] <= 1, "WRK fixed input is not a bit");
      fixed_masked[w] = fixed_bits[w] ^ masks[w];
    } else late_masks.push_back(masks[w]);
  }
  size_t ai = 0;
  for (const auto &g : program.gates) {
    using circuit::Op;
    if (g.op == Op::And) {
      const uint8_t ab = masks[g.in0] & masks[g.in1];
      prg.random_data_unaligned(&masks[g.out], 1); masks[g.out] &= 1;
      const auto out = test_dealer_detail::share_bit<nP>(masks[g.out], materials, prg);
      const auto product = test_dealer_detail::share_bit<nP>(ab, materials, prg);
      for (int p = 1; p <= nP; ++p) { materials[p].and_masks[ai] = out[p]; materials[p].products[ai] = product[p]; }
      ++ai;
    } else if (g.op == Op::Xor) masks[g.out] = masks[g.in0] ^ masks[g.in1];
    else if (g.op == Op::Not) masks[g.out] = masks[g.in0];
    else masks[g.out] = 0;
  }
  std::array<GcState<nP>, nP + 1> states;
  for (int p = 1; p <= nP; ++p)
    states[p] = prepare_local(program, materials[p], seed, fixed_inputs, fixed_masked, late_masks);
  auto &ev = states[1];
  ev.rows.resize(size_t(4) * (nP - 1) * na * row_bytes<nP>());
  for (int p = 2; p <= nP; ++p) {
    ev.fixed_input_labels[p].resize(fixed_inputs);
    for (uint32_t w = 0; w < fixed_inputs; ++w)
      ev.fixed_input_labels[p][w] = states[p].input_label0[w] ^
          (select_mask[fixed_masked[w]] & states[p].Delta);
    for (size_t a = 0; a < na; ++a)
      std::copy_n(states[p].rows.data() + 4 * a * row_bytes<nP>(),
                  4 * row_bytes<nP>(), ev.rows.data() + row_offset<nP>(a, p, 0));
    std::vector<uint8_t>().swap(states[p].rows); // rows already delivered offline
  }
  return states;
}

}  // namespace emp::wrk
