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
// Online, only the challenge-dependent late public inputs (the checked residues
// delta_H) drive one flight of selected labels to P1, who evaluates and decodes.
//
// The circuit itself is supplied as a compiled emp-tool BooleanProgram; this
// file is parameter- and mode-agnostic (both sign.h and sign_2round.h use it).
#ifndef MLDSA_WRK_PHASE2_H
#define MLDSA_WRK_PHASE2_H

#include "backend.h"

#include <emp-ag/wrk.h>
#include <emp-tool/ir/program.h>

#include <cstdint>
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
};

// ---- offline: sample, garble, deliver, open output masks --------------------
//
// fixed_shares are the fixed secret inputs as authenticated GMW shares, in the
// program's input order, exactly `fixed_inputs` of them. late_inputs =
// num_inputs - fixed_inputs are the public post-challenge ports.
template <int nP>
inline WrkOffline<nP> wrk_offline(Backend<nP>& bk, const emp::circuit::BooleanProgram& prog,
                                  const emp::wrk::ShareVec<nP>& fixed_shares,
                                  uint32_t fixed_inputs) {
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

  // 3. install fixed secret inputs (open B/M XOR lambda) and the late public
  //    port masks (opened to all) in ONE checked GMW open.
  emp::wrk::ShareVec<nP> bindings((size_t)n_in);
  for (uint32_t j = 0; j < fixed_inputs; ++j)
    bindings[j] = emp::wrk::xor_share(fixed_shares[j], plan.local.input_masks[j]);
  for (uint32_t j = fixed_inputs; j < n_in; ++j)
    bindings[j] = plan.local.input_masks[j];
  const std::vector<uint8_t> masked = bk.gmw().open(bindings, "WRK-fixed-and-late-masks");

  // 4. garble locally. Rows are produced by garblers; P1 keeps the bases.
  WrkOffline<nP> off;
  off.fixed_inputs = fixed_inputs;
  off.late_inputs = late;
  off.outputs = (uint32_t)prog.outputs.size();
  off.prog = &prog;
  off.gc = emp::wrk::prepare_local<nP>(prog, plan.local, bk.session_id(), fixed_inputs,
                                       std::span(masked).first(fixed_inputs),
                                       std::span(masked).subspan(fixed_inputs));

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
  return off;
}

// ---- online: one flight of selected labels, then local eval + decode --------
//
// late_bits are the public post-challenge input bits (the checked delta_H
// residues), in program input order after the fixed ports. Returns P1's decoded
// output bits; garblers get nullopt.
template <int nP>
inline std::optional<std::vector<uint8_t>>
wrk_online(Backend<nP>& bk, WrkOffline<nP>& off, const std::vector<uint8_t>& late_bits) {
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

} // namespace mldsa

#endif // MLDSA_WRK_PHASE2_H
