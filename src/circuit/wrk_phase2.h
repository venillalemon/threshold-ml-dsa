// src/circuit/wrk_phase2.h — the WRK17b online circuit (C_post), garbled OFFLINE.
//
// This realizes the paper's Phase-2/Phase-3 split with two post-challenge
// flights. Everything expensive is pre-challenge:
//   * GMW samples the WRK wire masks and the authenticated products
//     lambda_a*lambda_b for every AND (sample_wrk_masks + and_batch);
//   * the fixed secret inputs (B, M, and in slot mode the compaction outputs
//     and rho) are installed by opening their differences from fresh GC masks;
//   * every garbler garbles its four rows and ships them, its fixed-input
//     labels, AND its output-mask shares to the evaluator P1 — so the output
//     is decoded locally online with no acceptance-conditioned release.
// Online there is exactly one garbler->P1 flight after the challenge:
//   * two-round mode: the late ports are PUBLIC (delta_H); garblers send the
//     selected labels (wrk_online).
//   * slot mode: the late ports are the routed residues d^ = XOR of retained u
//     wires with PUBLIC coefficients (delta_H). P1 forms d^'s active labels and
//     masked bits itself by free XOR; the ports themselves carry fresh offline
//     masks, so every gate is garbled offline, and flight 2 only carries the
//     re-masking from the XORed u wires onto the ports: one 1+16n-byte block
//     per port per garbler, the same content as one garbled row without the
//     padding (wrk_online_routed).
//
// The circuit itself is supplied as a compiled emp-tool BooleanProgram; this
// file is parameter- and mode-agnostic (sign_slot.h and sign_2round.h use it).
#ifndef MLDSA_WRK_PHASE2_H
#define MLDSA_WRK_PHASE2_H

#include "backend.h"

#include <emp-ag/wrk.h>
#include <emp-tool/ir/program.h>

#include <cstdint>
#include <cstring>
#include <future>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace mldsa {

// Concurrent all-to-all: run send(peer) then recv(peer) for every peer, in
// parallel over the pool, draining before returning. Safe on the quiescent
// mesh between GMW/arithmetic phases (symmetric, fully joined).
template <int nP, class Send, class Recv>
inline void mesh_exchange(Backend<nP>& bk, Send&& send, Recv&& recv) {
  std::vector<std::future<void>> jobs;
  for (int p = 1; p <= nP; ++p)
    if (p != bk.party()) {
      jobs.push_back(bk.pool()->enqueue([&, p] { send(p); bk.io().flush(p); }));
      jobs.push_back(bk.pool()->enqueue([&, p] { recv(p); }));
    }
  for (auto& j : jobs)
    j.get();
}

// A WRK preprocessing plan plus the AND-input mask pairs whose products the
// caller closes with one GMW and_batch.
template <int nP> struct MaskPlan {
  emp::wrk::LocalPreprocessing<nP> local;
  emp::wrk::ShareVec<nP> left, right;
};

// Sample authenticated WRK wire masks through GMW, then walk the circuit to
// gather the actual AND-input mask pairs for one explicit product stage. The
// caller runs plan.local.products = gmw.and_batch(plan.left, plan.right).
// Only garblers sample label-zero blocks (their own PRG); P1 samples none.
template <int nP>
inline MaskPlan<nP> wrk_sample_masks(Backend<nP>& bk, const emp::circuit::BooleanProgram& prog,
                                     int party, emp::block Delta) {
  emp::circuit::validate_program(prog);
  const uint32_t na = emp::wrk::count_gc_ands(prog);
  const size_t total = (size_t)prog.num_inputs + na;
  auto sampled = bk.gmw().random(total);

  MaskPlan<nP> plan;
  plan.local.party = party;
  plan.local.Delta = Delta;
  plan.local.input_masks.assign(sampled.begin(), sampled.begin() + prog.num_inputs);
  plan.local.and_masks.assign(sampled.begin() + prog.num_inputs, sampled.end());
  plan.left.reserve(na);
  plan.right.reserve(na);

  emp::wrk::ShareVec<nP> wires((size_t)prog.num_wires);
  std::copy(plan.local.input_masks.begin(), plan.local.input_masks.end(), wires.begin());
  const auto zero = bk.gmw().public_bit(false);
  size_t idx = 0;
  for (const auto& g : prog.gates) {
    using emp::circuit::Op;
    switch (g.op) {
    case Op::Xor:
      wires[g.out] = emp::wrk::xor_share(wires[g.in0], wires[g.in1]);
      break;
    case Op::Not:
      wires[g.out] = wires[g.in0];
      break;
    case Op::Const0:
    case Op::Const1:
      wires[g.out] = zero;
      break;
    case Op::And:
      plan.left.push_back(wires[g.in0]);
      plan.right.push_back(wires[g.in1]);
      wires[g.out] = plan.local.and_masks[idx++];
      break;
    default:
      throw std::runtime_error("wrk_sample_masks: unsupported opcode");
    }
  }
  if (party != 1) {
    plan.local.input_labels0.resize(prog.num_inputs);
    plan.local.and_labels0.resize(na);
    emp::PRG rng;
    if (!plan.local.input_labels0.empty())
      rng.random_block(plan.local.input_labels0.data(), plan.local.input_labels0.size());
    if (!plan.local.and_labels0.empty())
      rng.random_block(plan.local.and_labels0.data(), plan.local.and_labels0.size());
  }
  return plan;
}

// Offline state for one prepared C_post execution.
template <int nP> struct WrkOffline {
  emp::wrk::GcState<nP> gc;             // P1: full table + bases; garbler: local rows
  uint32_t fixed_inputs = 0;            // number of pre-challenge secret input ports
  uint32_t late_inputs = 0;            // number of post-challenge public input ports
  uint32_t outputs = 0;
  std::vector<uint8_t> out_mask_opened; // P1 only: XOR of all parties' output masks
  const emp::circuit::BooleanProgram* prog = nullptr;
  // Late ports are PUBLIC (masks opened offline, wrk_online) unless
  // late_routed: then they are formed online by free XOR of retained fixed
  // wires and their placeholder masks are never opened (wrk_online_routed).
  bool late_routed = false;
  // Retained WRK preprocessing: everyone's wire-mask shares and, at garblers,
  // the zero labels. wrk_online_routed needs input_masks (u and d^ ports) and
  // input_labels0 (garblers) to build the re-masking blocks once the routing
  // is known.
  emp::wrk::LocalPreprocessing<nP> local;
};

// ---- offline: sample, garble, deliver, open output masks --------------------
//
// fixed_shares are the fixed secret inputs as authenticated GMW shares, in the
// program's input order, exactly `fixed_inputs` of them. late_inputs =
// num_inputs - fixed_inputs are the public post-challenge ports.
template <int nP>
inline WrkOffline<nP> wrk_offline(Backend<nP>& bk, const emp::circuit::BooleanProgram& prog,
                                  const emp::wrk::ShareVec<nP>& fixed_shares,
                                  uint32_t fixed_inputs, bool late_routed = false) {
  using emp::wrk::AuthShare;
  const int party = bk.party();
  const emp::block Delta = bk.delta();
  const uint32_t na = emp::wrk::count_gc_ands(prog);
  const uint32_t n_in = prog.num_inputs;
  emp::expecting(fixed_inputs <= n_in && fixed_shares.size() == fixed_inputs,
                 "wrk_offline: fixed input width mismatch");
  const uint32_t late = n_in - fixed_inputs;

  // 1. sample wire masks + gather AND-input mask pairs, 2. products via GMW.
  auto plan = wrk_sample_masks(bk, prog, party, Delta);
  plan.local.products = bk.gmw().and_batch(plan.left, plan.right);
  plan.left.clear();
  plan.right.clear();

  // 3. install fixed secret inputs (open B/M XOR lambda) and, for PUBLIC late
  //    ports, their masks (opened to all) -- ONE checked GMW open. ROUTED late
  //    ports get their real masks online (XOR of retained fixed-wire masks);
  //    the placeholders sampled here are never opened.
  const uint32_t n_open = late_routed ? fixed_inputs : n_in;
  emp::wrk::ShareVec<nP> bindings((size_t)n_open);
  for (uint32_t j = 0; j < fixed_inputs; ++j)
    bindings[j] = emp::wrk::xor_share(fixed_shares[j], plan.local.input_masks[j]);
  if (!late_routed)
    for (uint32_t j = fixed_inputs; j < n_in; ++j)
      bindings[j] = plan.local.input_masks[j];
  const std::vector<uint8_t> masked = bk.gmw().open(bindings, "WRK-fixed-and-late-masks");
  std::vector<uint8_t> late_masks;
  if (late_routed)
    late_masks.assign((size_t)late, 0); // placeholder; online masked bits are set directly
  else
    late_masks.assign(masked.begin() + fixed_inputs, masked.end());

  // 4. garble locally. Rows are produced by garblers; P1 keeps the bases.
  WrkOffline<nP> off;
  off.fixed_inputs = fixed_inputs;
  off.late_inputs = late;
  off.outputs = (uint32_t)prog.outputs.size();
  off.prog = &prog;
  off.late_routed = late_routed;
  off.gc = emp::wrk::prepare_local<nP>(prog, plan.local, bk.session_id(), fixed_inputs,
                                       std::span(masked).first(fixed_inputs),
                                       std::span(late_masks));

  // 5. deliver rows + fixed-input labels to P1 (garblers send, P1 assembles).
  const size_t rb = emp::wrk::row_bytes<nP>();
  if (party == 1) {
    off.gc.rows.resize((size_t)4 * (nP - 1) * na * rb);
    std::vector<std::future<void>> jobs;
    for (int peer = 2; peer <= nP; ++peer)
      jobs.push_back(bk.pool()->enqueue([&, peer] {
        std::vector<uint8_t> rows((size_t)4 * na * rb);
        if (na)
          bk.io().recv_data(peer, rows.data(), rows.size());
        for (uint32_t a = 0; a < na; ++a)
          std::copy_n(rows.data() + (size_t)4 * a * rb, 4 * rb,
                      off.gc.rows.data() + emp::wrk::row_offset<nP>(a, peer, 0));
        off.gc.fixed_input_labels[peer].resize(fixed_inputs);
        if (fixed_inputs)
          bk.io().recv_data(peer, off.gc.fixed_input_labels[peer].data(),
                            (size_t)fixed_inputs * sizeof(emp::block));
      }));
    for (auto& j : jobs)
      j.get();
  } else {
    emp::BlockVec sel((size_t)fixed_inputs);
    for (uint32_t j = 0; j < fixed_inputs; ++j)
      sel[j] = off.gc.input_label0[j] ^ (emp::select_mask[off.gc.input_masked[j]] & Delta);
    if (na)
      bk.io().send_data(1, off.gc.rows.data(), off.gc.rows.size());
    if (fixed_inputs)
      bk.io().send_data(1, sel.data(), (size_t)fixed_inputs * sizeof(emp::block));
    bk.io().flush(1);
    std::vector<uint8_t>().swap(off.gc.rows); // delivered
  }

  // 6. open the output masks to P1 (offline): garblers send bit + MAC, P1
  //    verifies each and forms the XOR mask, so online decode is local.
  const uint32_t no = off.outputs;
  if (party != 1) {
    std::vector<uint8_t> bits((size_t)no);
    emp::BlockVec macs((size_t)no);
    for (uint32_t i = 0; i < no; ++i) {
      bits[i] = off.gc.output_masks[i].bit;
      macs[i] = off.gc.output_masks[i].mac(emp::wrk::peer_slot(party, 1));
    }
    bk.io().send_data(1, bits.data(), bits.size());
    if (no)
      bk.io().send_data(1, macs.data(), (size_t)no * sizeof(emp::block));
    bk.io().flush(1);
  } else {
    off.out_mask_opened.assign((size_t)no, 0);
    for (uint32_t i = 0; i < no; ++i)
      off.out_mask_opened[i] = off.gc.output_masks[i].bit;
    std::vector<std::future<void>> jobs;
    std::vector<std::vector<uint8_t>> pbits((size_t)nP + 1);
    std::vector<emp::BlockVec> pmacs((size_t)nP + 1);
    for (int peer = 2; peer <= nP; ++peer) {
      pbits[peer].resize((size_t)no);
      pmacs[peer].resize((size_t)no);
      jobs.push_back(bk.pool()->enqueue([&, peer] {
        bk.io().recv_data(peer, pbits[peer].data(), (size_t)no);
        if (no)
          bk.io().recv_data(peer, pmacs[peer].data(), (size_t)no * sizeof(emp::block));
      }));
    }
    for (auto& j : jobs)
      j.get();
    for (int peer = 2; peer <= nP; ++peer)
      for (uint32_t i = 0; i < no; ++i) {
        emp::expecting(pbits[peer][i] <= 1, "wrk_offline: noncanonical output-mask bit");
        const emp::block expect =
            off.gc.output_masks[i].key(emp::wrk::peer_slot(1, peer)) ^
            (emp::select_mask[pbits[peer][i]] & bk.delta());
        emp::expecting(emp::cmpBlock(&expect, &pmacs[peer][i], 1),
                       "wrk_offline: output-mask MAC verification failed");
        off.out_mask_opened[i] ^= pbits[peer][i];
      }
  }
  // Retain the preprocessing the routed-online path needs (products included).
  off.local = std::move(plan.local);
  return off;
}

// ---- online: one flight of selected labels, then local eval + decode --------
//
// late_bits are the (public) masked bits of the late ports, in program input
// order after the fixed ports: for PUBLIC late ports the bit itself (its mask
// was opened offline), for SECRET late ports the opened value ^ mask. Returns
// P1's decoded output bits; garblers get nullopt.
template <int nP>
inline std::optional<std::vector<uint8_t>>
wrk_online_bits(Backend<nP>& bk, WrkOffline<nP>& off, const std::vector<uint8_t>& late_bits) {
  const int party = bk.party();
  emp::expecting(late_bits.size() == off.late_inputs, "wrk_online: late input width mismatch");
  const uint32_t n_in = off.fixed_inputs + off.late_inputs;

  if (party != 1) {
    // garbler: select the labels for the public late ports and send to P1.
    emp::BlockVec sel = emp::wrk::selected_late_labels(off.gc, std::span(late_bits));
    if (!sel.empty())
      bk.io().send_data(1, sel.data(), sel.size() * sizeof(emp::block));
    bk.io().flush(1);
    return std::nullopt;
  }

  // evaluator: bind public ports, receive each garbler's late labels, evaluate.
  emp::wrk::bind_public_inputs(off.gc, std::span(late_bits));
  std::array<emp::BlockVec, nP + 1> labels;
  std::vector<std::future<void>> jobs;
  for (int peer = 2; peer <= nP; ++peer) {
    labels[peer] = off.gc.fixed_input_labels[peer];
    labels[peer].resize((size_t)n_in);
    jobs.push_back(bk.pool()->enqueue([&, peer] {
      if (off.late_inputs)
        bk.io().recv_data(peer, labels[peer].data() + off.fixed_inputs,
                          (size_t)off.late_inputs * sizeof(emp::block));
    }));
  }
  for (auto& j : jobs)
    j.get();

  std::vector<uint8_t> masked = emp::wrk::evaluate(off.gc, *off.prog, labels, bk.pool());
  emp::expecting(masked.size() == off.out_mask_opened.size(), "wrk_online: output width mismatch");
  std::vector<uint8_t> out((size_t)off.outputs);
  for (uint32_t i = 0; i < off.outputs; ++i)
    out[i] = (uint8_t)(masked[i] ^ off.out_mask_opened[i]);
  return out;
}

// PUBLIC late ports (two-round mode): the late values are public bits.
template <int nP>
inline std::optional<std::vector<uint8_t>>
wrk_online(Backend<nP>& bk, WrkOffline<nP>& off, const std::vector<uint8_t>& late_bits) {
  emp::expecting(!off.late_routed, "wrk_online: this circuit has ROUTED late ports");
  return wrk_online_bits<nP>(bk, off, late_bits);
}

// ---- online, ROUTED late ports (slot mode) ----------------------------------
//
// Port layout the caller guarantees: u[j][i] is fixed port u_base + j*M + i
// (installed offline, never read by a gate) and d^[j][l] is late port
// dh_base + j*K + l, read by the offline-garbled gates under its own fresh
// mask lambda_f (sampled offline like any input port). S[l] lists the
// coefficients i with delta_H,i[l] = 1 (public after flight 1).
//
// d^[j][l] = XOR_{i in S[l]} u[j][i]. By free XOR P1 holds, for every garbler
// i, the active label L_i(d) = L_i,0(d) ^ Lambda_d * Delta_i and the masked
// bit Lambda_d = x ^ lambda_d, where lambda_d = XOR of the u masks (shared).
// What the offline-garbled gates need is the same x under the port's mask:
//   Lambda_f = Lambda_d ^ delta,               delta = lambda_d ^ lambda_f,
//   L_i(f)   = L_i(d) ^ [L_i,0(f) ^ L_i,0(d) ^ delta * Delta_i].
// delta is only ever shared, so P1 assembles delta*Delta_i from the MACs the
// way evaluate() assembles a row: each garbler i sends, per port,
//   [bit of its delta share][its MACs on that share, one per other party]
//   [L_i,0(f) ^ L_i,0(d) ^ bit*Delta_i ^ XOR of its keys on the others' shares]
// = 1 + 16*nP bytes, the plaintext of one WRK row (no padding: there is a
// single block per port and P1 is entitled to it). P1 checks every bit
// against its own key, so a lying garbler is caught here; a wrong label is
// caught by the next row's MAC check. One garbler->P1 flight.
template <int nP>
inline std::optional<std::vector<uint8_t>>
wrk_online_routed(Backend<nP>& bk, WrkOffline<nP>& off, uint32_t u_base, uint32_t M,
                  uint32_t dh_base, uint32_t K, uint32_t NS,
                  const std::vector<std::vector<int>>& S) {
  using emp::wrk::AuthShare;
  using emp::wrk::peer_slot;
  using emp::wrk::xor_share;
  const int party = bk.party();
  const emp::block Delta = bk.delta();
  const uint32_t NB = NS * K; // routed ports
  emp::expecting(off.late_routed && off.late_inputs == NB && dh_base == off.fixed_inputs &&
                     S.size() == K,
                 "wrk_online_routed: port layout mismatch");
  constexpr size_t blk = 1 + size_t(16) * nP; // one re-masking block

  // Everyone: own share of delta = lambda_d ^ lambda_f per port.
  emp::wrk::ShareVec<nP> dsh((size_t)NB);
  for (uint32_t j = 0; j < NS; ++j)
    for (uint32_t l = 0; l < K; ++l) {
      const uint32_t w = j * K + l;
      AuthShare<nP> a = off.local.input_masks[dh_base + w];
      for (int i : S[l])
        a = xor_share(a, off.local.input_masks[u_base + j * M + (uint32_t)i]);
      dsh[w] = a;
    }

  if (party != 1) {
    std::vector<uint8_t> out((size_t)NB * blk);
    for (uint32_t j = 0; j < NS; ++j)
      for (uint32_t l = 0; l < K; ++l) {
        const uint32_t w = j * K + l;
        emp::block corr = off.local.input_labels0[dh_base + w]; // L_0(f)
        for (int i : S[l])
          corr = corr ^ off.local.input_labels0[u_base + j * M + (uint32_t)i]; // ^ L_0(d)
        const AuthShare<nP>& sh = dsh[w];
        uint8_t* b = out.data() + (size_t)w * blk;
        b[0] = sh.bit;
        corr = corr ^ (emp::select_mask[sh.bit] & Delta);
        for (int peer = 1; peer <= nP; ++peer)
          if (peer != party) {
            const int slot = peer_slot(party, peer);
            emp::wrk::detail::put_block(b + 1 + 16 * slot, sh.mac(slot));
            corr = corr ^ sh.key(slot);
          }
        emp::wrk::detail::put_block(b + 1 + 16 * (nP - 1), corr);
      }
    bk.io().send_data(1, out.data(), out.size());
    bk.io().flush(1);
    return std::nullopt;
  }

  // P1: receive every garbler's blocks.
  std::array<std::vector<uint8_t>, nP + 1> in;
  std::vector<std::future<void>> jobs;
  for (int peer = 2; peer <= nP; ++peer) {
    in[peer].resize((size_t)NB * blk);
    jobs.push_back(bk.pool()->enqueue(
        [&, peer] { bk.io().recv_data(peer, in[peer].data(), in[peer].size()); }));
  }
  for (auto& jb : jobs)
    jb.get();

  // P1: d^ by free XOR, then re-mask onto the ports.
  const uint32_t n_in = off.fixed_inputs + off.late_inputs;
  std::array<emp::BlockVec, nP + 1> labels;
  for (int peer = 2; peer <= nP; ++peer) {
    labels[peer] = off.gc.fixed_input_labels[peer];
    labels[peer].resize((size_t)n_in);
  }
  for (uint32_t j = 0; j < NS; ++j)
    for (uint32_t l = 0; l < K; ++l) {
      const uint32_t w = j * K + l, port = dh_base + w;
      uint8_t lam = 0;
      std::array<emp::block, nP + 1> lab{};
      for (int i : S[l]) {
        const uint32_t up = u_base + j * M + (uint32_t)i;
        lam ^= off.gc.input_masked[up];
        for (int peer = 2; peer <= nP; ++peer)
          lab[peer] = lab[peer] ^ off.gc.fixed_input_labels[peer][up];
      }
      // delta: own share bit + every garbler's, each checked under P1's key.
      const AuthShare<nP>& mine = dsh[w];
      uint8_t delta = mine.bit;
      for (int peer = 2; peer <= nP; ++peer) {
        const uint8_t* b = in[peer].data() + (size_t)w * blk;
        emp::expecting(b[0] <= 1, "wrk_online_routed: noncanonical re-mask bit");
        const emp::block expect = mine.key(peer_slot(1, peer)) ^ (emp::select_mask[b[0]] & Delta);
        const emp::block got = emp::wrk::detail::get_block(b + 1 + 16 * peer_slot(peer, 1));
        emp::expecting(emp::cmpBlock(&expect, &got, 1),
                       "wrk_online_routed: re-mask MAC verification failed");
        delta ^= b[0];
      }
      off.gc.input_masked[port] = (uint8_t)(lam ^ delta);
      // labels: L(d) ^ corr_dest ^ (all MACs on the others' shares under Delta_dest)
      for (int dest = 2; dest <= nP; ++dest) {
        const uint8_t* bd = in[dest].data() + (size_t)w * blk;
        emp::block v = emp::wrk::detail::get_block(bd + 1 + 16 * (nP - 1));
        v = v ^ mine.mac(peer_slot(1, dest));
        for (int src = 2; src <= nP; ++src)
          if (src != dest)
            v = v ^ emp::wrk::detail::get_block(in[src].data() + (size_t)w * blk + 1 +
                                                16 * peer_slot(src, dest));
        labels[dest][port] = lab[dest] ^ v;
      }
    }

  std::vector<uint8_t> masked = emp::wrk::evaluate(off.gc, *off.prog, labels, bk.pool());
  emp::expecting(masked.size() == off.out_mask_opened.size(), "wrk_online_routed: output width");
  std::vector<uint8_t> out((size_t)off.outputs);
  for (uint32_t i = 0; i < off.outputs; ++i)
    out[i] = (uint8_t)(masked[i] ^ off.out_mask_opened[i]);
  return out;
}

} // namespace mldsa

#endif // MLDSA_WRK_PHASE2_H
