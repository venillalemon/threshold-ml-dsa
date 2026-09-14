#ifndef EMP_AG_BACKEND_PREPROC_BUFFER_OT_H__
#define EMP_AG_BACKEND_PREPROC_BUFFER_OT_H__

// BufferOT: realize a base OT from PRECOMPUTED random OTs via Beaver's OT
// precomputation, so each unordered party pair needs only one public-key base
// OT for its two COT directions. The precomputed OTs come from the canonical
// primary direction (128 COTs, CCR-hashed to break the Delta correlation).
// Implements emp-ot's OT interface and uses every OTExt's common fourth
// `base_ot` constructor argument. Party-specific: the extension-SENDER's base
// OT is a Beaver RECEIVER (holds (choice, m_choice)); the
// extension-RECEIVER's is a Beaver SENDER (holds (m0, m1)).
//
// Beaver correction (precomputed sender (m0,m1), receiver (c, m_c); realize a
// chosen OT with sender inputs (K0,K1), receiver choice d):
//   receiver: e = d ^ c ; send e
//   sender  : Y_j = K_j ^ m_{j^e}  (j in {0,1}) ; send Y0,Y1
//   receiver: out = Y_d ^ m_c      [ = K_d, since d^e = c ]
//
// Correctness + malicious base-OT soundness validated standalone: SS#2 seeded
// this way produces correct, checkable COTs.

#include "emp-ot/emp-ot.h"
#include <cstdint>
#include <utility>
#include <vector>

namespace emp::ag {

class BufferOT : public emp::OT {
 public:
  BufferOT(emp::IOChannel* io, std::vector<uint8_t> choice, std::vector<block> m_choice)
      : io_(io), c_(std::move(choice)), mc_(std::move(m_choice)) {
    expecting(io_ != nullptr, "BufferOT: IO channel must not be null");
    expecting(c_.size() == mc_.size(),
              "BufferOT: receiver choice/message buffers must have equal length");
  }
  BufferOT(emp::IOChannel* io, std::vector<block> m0, std::vector<block> m1)
      : io_(io), m0_(std::move(m0)), m1_(std::move(m1)) {
    expecting(io_ != nullptr, "BufferOT: IO channel must not be null");
    expecting(m0_.size() == m1_.size(),
              "BufferOT: sender message buffers must have equal length");
  }

  bool is_malicious_secure() const override { return true; }

  void send(const block* K0, const block* K1, int64_t n) override {
    expecting(n >= 0, "BufferOT::send: negative OT count");
    expecting(n == 0 || (K0 != nullptr && K1 != nullptr),
              "BufferOT::send: null input for nonzero count");
    expecting(pos_ <= m0_.size() && (size_t)n <= m0_.size() - pos_,
              "BufferOT::send: precomputed OT buffer exhausted");
    std::vector<uint8_t> e((size_t)n);
    io_->recv_data(e.data(), (size_t)n);
    std::vector<block> Y(2 * (size_t)n);
    for (int64_t i = 0; i < n; ++i) {
      const size_t j = pos_ + (size_t)i;
      const block m_e  = e[(size_t)i] ? m1_[j] : m0_[j];
      const block m_ne = e[(size_t)i] ? m0_[j] : m1_[j];
      Y[2 * (size_t)i]     = K0[i] ^ m_e;
      Y[2 * (size_t)i + 1] = K1[i] ^ m_ne;
    }
    io_->send_data(Y.data(), 2 * (size_t)n * sizeof(block));
    io_->flush();
    pos_ += (size_t)n;
  }

  void recv(block* out, const bool* b, int64_t n) override {
    expecting(n >= 0, "BufferOT::recv: negative OT count");
    expecting(n == 0 || (out != nullptr && b != nullptr),
              "BufferOT::recv: null input/output for nonzero count");
    expecting(pos_ <= c_.size() && (size_t)n <= c_.size() - pos_,
              "BufferOT::recv: precomputed OT buffer exhausted");
    std::vector<uint8_t> e((size_t)n);
    for (int64_t i = 0; i < n; ++i)
      e[(size_t)i] = (uint8_t)(b[i] ^ (c_[pos_ + (size_t)i] & 1));
    io_->send_data(e.data(), (size_t)n);
    io_->flush();
    std::vector<block> Y(2 * (size_t)n);
    io_->recv_data(Y.data(), 2 * (size_t)n * sizeof(block));
    for (int64_t i = 0; i < n; ++i)
      out[i] = (b[i] ? Y[2 * (size_t)i + 1] : Y[2 * (size_t)i]) ^
               mc_[pos_ + (size_t)i];
    pos_ += (size_t)n;
  }

 private:
  emp::IOChannel* io_;
  std::vector<uint8_t> c_;
  std::vector<block>   mc_;
  std::vector<block>   m0_, m1_;
  size_t pos_ = 0;
};

}  // namespace emp::ag
#endif  // EMP_AG_BACKEND_PREPROC_BUFFER_OT_H__
