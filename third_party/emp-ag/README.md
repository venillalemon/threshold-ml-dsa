# Local classical WRK backend

This is the repository-local authenticated-garbling backend for
[`two-round-online`](../two-round-online/). It implements the four-row
construction in WRK17b, Figures 2–3 of
[Global-Scale Secure Multiparty Computation](https://eprint.iacr.org/2017/189).
It is not a full fork of the sibling KRRW library. That checkout is unchanged.

The interface is in [`emp-ag/wrk.h`](emp-ag/wrk.h), namespace `emp::wrk`.
The CMake target is `emp-ag::wrk`, with dependencies on `emp-tool` and
OpenSSL. No OT initialization or networking occurs inside the garbler or
evaluator. The separately copied `backend/netmp.h` is a test transport;
its provenance and license are recorded in
[`DEPENDENCIES.md`](../two-round-online/DEPENDENCIES.md).

The separate [`emp-ag/gmw.h`](emp-ag/gmw.h) adapter provides real
OT-based authenticated Boolean GMW through target `emp-ag::gmw`, which
also requires `emp-ot`. It is used by the distributed preprocessing
builder; it does not add online interaction inside the WRK evaluator.
Its default COT extension is Ferret; CMake option
`EMP_AG_OT_BACKEND=SoftSpoken` retains the previous backend for comparison,
and `IKNP` is also available. Ferret adds an LPN assumption and still uses
SoftSpoken internally for bootstrap. All malicious-COT, authenticated
opening, and checked-triple checks remain enabled. This changes
preprocessing communication, not the WRK construction or online schedule.

`EMP_AG_TRIPLE_BACKEND=FKOS` optionally replaces the default `HalfGate`
random-triple provider with FKOS cut-and-choose/double bucketing. It uses
the old authenticated-bit/hash flow on the same COT streams and Delta,
with no extra edaBit assumptions. Its default profile uses large batches
and 40-bit bucket security. See [FKOS.md](FKOS.md) for checks, the hash-reuse
proof limitation, parameters, and tests.

## Preprocessing contract

Each party has its own private 128-bit verifier key `Delta_i`. A garbler
uses that same `Delta_i` as its free-XOR offset. The keys are not a single
offset disclosed to all parties. Authentication stores the bit explicitly:

```
M_j[x_i] = K_j[x_i] XOR x_i * Delta_j.
```

No MAC/key bits are reserved for point-and-permute. `prepare_local` takes
only one party's `LocalPreprocessing`: authenticated input masks, fresh
AND-output masks, authenticated products of the actual AND-input masks,
and that garbler's independent zero-labels. These correlations must come
from secure preprocessing under the same verifier keys. Independent
Beaver triples cannot simply be substituted for the mask products.

Before using the returned state, preprocessing must also:

- Agree on the exact circuit and fresh hash seed.
- Install fixed secret inputs by opening their differences from fresh GC
  masks, and deliver only the corresponding active labels to the evaluator.
- Open masks for the late **public** ports; never open secret-port masks.
- Deliver every garbler's four rows per AND to the evaluator.
- Retain authenticated output-mask shares for the authorized output policy.

`make_test_material` supplies this contract using a centralized trusted
dealer for tests. It is not malicious-secure distributed preprocessing.
Fresh per-attempt state is mandatory. `prepare_local` itself cannot certify
the consistency of externally supplied correlations.

`two-round-online/mldsa_preprocess` realizes this preparation through
party-separated GMW and table delivery, assuming free authenticated
eDaBits and key-generation bits. Its free test provider supplies no
triples, boundary values, wire masks, mask products, or garbled tables.
The native GMW representation fixes two low offset bits, so that path has
126-bit Boolean key entropy despite full 128-bit blocks. Checked triple
generation uses statistical parameter 80 per batch; a union bound over
`B` relevant batches is `B*2^-80`, not a whole-session guarantee of 80
bits. See [PREPROCESSING.md](../two-round-online/PREPROCESSING.md).

## Row construction and local evaluation

For an AND gate with wire masks `a,b,g`, form the authenticated affine row

```
r_uv = a*b XOR g XOR u*b XOR v*a XOR (u AND v).
```

The public last term is injected into party 1's share with the corresponding
verifier-key adjustment. Garbler `i` encrypts its row share, its MACs to all
other parties, and

```
L_g,0^i XOR r_uv^i * Delta_i XOR XOR_{j != i} K_i[r_uv^j].
```

The pad is SHAKE256 on the ordered pair of active input labels, with domain
separation for the execution seed, gate, row, garbler, and output length.
This instantiates the paper's hash in the random-oracle model; it is not
an unaudited substitution of a half-gate hash.

The evaluator decrypts one row per garbler, verifies every row share's MAC
under its private key, and then reconstructs the masked output bit and
all active output labels. XOR and NOT gates propagate locally. Evaluation
returns **masked** output bits or throws; authenticated output decoding is
handled by the caller. There is no KRRW witness broadcast, joint coin, or
interactive zero check.

For seven parties, each row is `1 + 7*16 = 113` bytes. The raw table is
`4*6*113 = 2,712` bytes per AND, or **147,182,952 bytes** for the unchanged
54,271-AND MLDSA44 Phase-2 circuit. The evaluator additionally stores three
authenticated affine bases per AND: **31,422,909 serialized bytes**. These
are pre-challenge costs, not online communication. This implementation
does not apply row reduction or ciphertext aggregation.

## Validation and scope

Build and test via `two-round-online/CMakeLists.txt`. Tests cover two and
seven parties, every small-circuit input assignment, all-public/all-fixed
and mixed input boundaries, constants, free gates, direct input outputs,
corrupted row MACs and labels, and the production circuit. Network tests
exercise the two-flight application schedule separately.

The two-flight claim is for **evaluator-local output with abort**, after
preprocessing and a common challenge, over assumed authenticated private
channels. It does not include reliable broadcast, fairness, simultaneous
common abort, or output dissemination. See
[`SECURITY.md`](../two-round-online/SECURITY.md) for these boundaries.
This research implementation and its tests are not a security audit.
