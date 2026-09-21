# threshold-ml-dsa

Threshold ML-DSA signing (FIPS 204) with a mixed boolean backend: the
pre-challenge circuit C_pre runs under authenticated GMW (TinyOT triples), the
post-challenge circuit C_post is a WRK17 four-row garbled circuit that is
garbled and delivered offline. Two modes, same offline/online structure, only
the circuits differ:

- **two-round** (`SLOT=0`): Π_TwoRound of the paper.
- **slot** (`SLOT=1`): Π_Slot; d^ is routed by free XOR with public δ_H, the
  d^ ports are re-masked online with one 1+16n-byte block per port per garbler.

Both modes need exactly two online flights after the message.

## Layout

```
CMakeLists.txt, Makefile      make n=<parties> p=<44|65|87> SLOT=<0|1> [TEST=1] [OT=Ferret]
main.cpp                      KeyGen + one Sign; prints gate counts and communication
runN                          local n-party launcher
src/infra/backend.h           NetIOMP mesh + GMW engine + per-party Δ
src/infra/bdoz.h              BDOZ shares over F_q and Z_{2^k}, checked one-flight opens
src/infra/dealer.h            demo dealer: keys, edaBits, Boolean shares (cheats, see below)
src/circuit/a2b.h, decompose.h, circuit.h
                              circuit pieces (A2B, Decompose, adders, MUX)
src/circuit/wrk_phase2.h      C_post glue: GMW masks -> WRK garbling -> delivery -> online
src/protocol/keygen.h         KeyGen
src/protocol/sign_2round.h    two-round Sign (C_pre + C_post + flights)
src/protocol/sign_slot.h      slot Sign
tests/                        backend, GMW and WRK bridge tests (MLDSA_TESTS=ON)
third_party/emp-ag/           vendored WRK/GMW backend (emp-ag::gmw); needs emp-tool, emp-ot
docs/                         Decompose correctness proof
```

Demo cheats (not part of the protocol cost): the dealer samples s, e, t, all
edaBits and their Boolean shares in the clear and hands out shares; every
party's Δ, session id and seed are deterministic. Everything after that is the
real protocol.

## Results

n = 7 parties, localhost, ssp = 80, σ = 128, SoftSpoken OT. Gates are AND
gates. "total" is the sum of bytes sent by all parties; "P1" is bytes sent
plus received by the evaluator.

|                        | ML-DSA-44 |        | ML-DSA-65 |        | ML-DSA-87 |        |
|------------------------|----------:|-------:|----------:|-------:|----------:|-------:|
|                        | two-round |   slot | two-round |   slot | two-round |   slot |
| offline gates (C_pre)  |   263,168 | 786,006 |  310,528 | 1,112,472 | 412,160 | 1,476,183 |
| online gates (C_post)  |    54,271 | 24,127 |    78,079 | 30,079 |   104,447 | 41,791 |
| offline comm, total (MB) | 1,747  |  4,145 |     2,175 |  5,832 |     2,891 |  7,757 |
| offline comm, P1 (MB)  |       608 |  1,243 |       778 |  1,741 |     1,036 |  2,320 |
| online comm, total (MB) |    3.244 |  2.037 |     4.745 |  2.690 |     6.081 |  3.455 |
| online comm, P1 (MB)   |     2.191 |  0.861 |     3.287 |  1.079 |     4.107 |  1.266 |

Online communication, two-round, ML-DSA-44: flight 1 is 42 × 35,104 B
(δ_H as 137-bit ring shares plus a digest, every ordered pair), flight 2 is
6 × 18,432 × 16 B (labels). Slot: flight 1 also carries V̄, flight 2 is
6 × 576 × 113 B (re-masking blocks).

With `OT=Ferret` offline communication drops by about 28% and offline time
grows about 3.5×.

## Build and run

```
make n=7 p=44 SLOT=0          # two-round, ML-DSA-44
make n=7 p=44 SLOT=1          # slot
make n=3 p=44 SLOT=1 TEST=1   # opens secrets, checks circuit against plaintext, verifies the signature
```

Requires emp-tool and emp-ot installed under `/usr/local`.
