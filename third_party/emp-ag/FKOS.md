# Large-batch cut-and-choose triple provider

`EMP_AG_TRIPLE_BACKEND=FKOS` selects the historical **CutChoose flow** in
[`fkos_triples.h`](emp-ag/backend/preproc/fkos_triples.h):

1. Use the existing COT/aShare provider to authenticate random `a,b,r`.
2. Hash `a`'s COT correlations into raw product shares `z`.
3. Open `z XOR r` and obtain authenticated `c=r XOR (z XOR r)`.
4. Check correctness with a random subset and bucket sacrifice.
5. Combine another set of buckets to remove leakage.

There is **one existing COT stream per directed peer pair and one existing
private Delta per party**. The hashing step reuses the authenticated `a`
correlations, as in the old implementation. No second stream, auxiliary
correlation batch, fresh Delta, or extra edaBit assumption is introduced.
`HalfGate` remains the default provider. MLDSA circuits, WRK rows, online
files, and the two online flights are unchanged.

Historical code references: `Documents/ref` commit `9a5e637` and
`Documents/emp-ag2pc` revision
`6096cd2^:emp-ag2pc/triple_pool_cutchoose.h`.

## Existing checks and completed checks

Candidate authentication calls the existing `AuthSharePool::draw` and
`maybe_flush_cot_check`. The native Delta-consistency check and aShare
choice/MAC-consistency check run as before. All COT checks complete before
hashing or opening masked products. There is no new cross-stream check.

Each COT-generation chunk appends 256 discarded authenticated random
bits. These mask the native aShare `check2`'s 128-bit linear disclosure:
the random binary 128-by-256 coefficient matrix has rank-deficiency
probability less than `2^-128`. The padding is never a candidate, a hash
input, or an output. This reuses the existing check rather than duplicating
or replacing its protocol. At two parties the native `check2` is vacuous;
the same padding layout is retained for simplicity.

The old prototype's missing opening checks are supplied by the existing
`checked_open_bits`: cut, sacrifice differences, sacrifice residuals, and
combination differences are all authenticated. MACs are checked as
domain-separated vector digests, and parties compare their opening views
and check status. Product-mask openings use echo broadcast.

Correctness and leakage buckets use independently sampled **uniform
permutations**, with unbiased Fisher-Yates shuffling. The correctness
proof is not extended to the old prototype's cyclic-shift shortcut.
Permutation coins are sampled after all candidates are generated. A
random-subset check rejects the all-bad-triples attack that pairwise
sacrifice alone cannot detect.

Candidate generation uses small chunks of 16,384 for bounded transient
COT buffers, but these are **not independent bucket batches**. All
candidates participate in one large global cut/sacrifice/combine batch.
Candidates are retained until the permutation is known; this uses more
RAM than the old 8,192-output batching. Surviving heads are compacted
in-place to avoid an additional full candidate copy.
The stream stays open across the authentication chunks; its malicious-COT
checks run once after the full population is authenticated, before the
hash/product pass. Chunk boundaries do not repeatedly finalize Ferret.

## Operations

For sender `j` and receiver `i`, the `a` authentication relation is
`M_j[a_i]=K_j[a_i] XOR a_i*Delta_j`. The sender sends one packed bit

```
H(K_j[a_i]) XOR H(K_j[a_i] XOR Delta_j) XOR b_j.
```

The receiver combines it with `H(M_j[a_i])` and `a_i`. Summing these
pairwise products and the local `a_i*b_i` terms gives shares of `a*b` in
an honest execution. Hashing uses emp-tool Hash (SHA-256 in the tested
build), separated by session, batch, ordered pair, and candidate index.

After authenticating `c`, correctness sacrifice against `(f,g,h)` opens
`rho=a XOR f`, `sigma=b XOR g`, and checks the authenticated zero value

```
c XOR h XOR rho*g XOR sigma*f XOR rho*sigma.
```

Leakage combination retains `b_head`, XORs all `a` values, and computes
`c_head XOR sum_other(c_other XOR (b_head XOR b_other)*a_other)`, using
authenticated openings of each `b_head XOR b_other`. Resulting random
triples enter the unchanged GMW Beaver bank. GMW still computes the actual
WRK mask products; random triples are not substituted for those products.

## Size-3 buckets at 40 bits

The current MLDSA44 preprocessing job consumes **322,559 triples**:
251,904 Phase-1 ANDs, 16,384 key-BtoA ANDs, and 54,271 WRK mask-product
ANDs. This fits one batch with correctness size `T=3`, combination size
`B=3`, and a random cut of 60 candidates. It generates **2,903,091 raw
candidates**, not 40 separate populations with larger buckets.

The parameter selector retains the finite-population factors. With
`ell` outputs and `t=B*ell` surviving heads, the bucket-analysis bounds are

```
correctness <= t / binomial(T*t,T),
leakage <= max_g 2^-g * ell*binomial(g,B) / binomial(B*ell,B).
```

The second maximum occurs at `g=2B-1` or `2B`. The first uses the
random-subset check with at least `T` cut items and the same correctness
bucketing argument as [FKOS15, Appendix F.3](https://eprint.iacr.org/2015/901.pdf).
The second is a first-moment bound for uniformly permuted leakage
buckets, using the leaky-triple survival bound `2^-g`. These are tighter
than replacing binomial coefficients by powers of the batch size.

At `ell=322559,T=B=3`, the bounds are approximately `2.37316e-13` and
`6.67453e-13`; their sum is `9.04769e-13 < 2^-40`. Exactly 250,000
outputs does not meet this combined 40-bit bound with both sizes equal
to 3. The selector increases a bucket size for such smaller batches.
This is a bucket-error bound, separately from authentication, COT, hash,
and padding-rank errors. Multiple bucket batches require a union bound.

For `n` candidates in one generation chunk, the explicit draw consumes
`3*n+256` authenticated bits per directed pair. Counters exclude the
underlying COT bootstrap/check internals; network byte counters include
those internal costs. `fkos_breakdown` is a subdivision of the outer
`gmw_triple_generation` stage, not an additional cost.

## Proof scope

This deliberately follows the historical **lean hash-reuse** construction,
not FKOS Figure 15's doubled independent correlations. FKOS footnote 4
says its UC proof does not cover reuse of the same outputs for both
authentication and hashing. Consequently, FKOS's theorem must not be
cited as a complete proof of this implementation. The concrete checks
and bucket bounds do not close that separate composition/proof question.
The provider remains an optional research implementation pending review;
functional and fault tests are not a cryptographic audit. The native
profile has 126 bits of Boolean key entropy, not 128.

## Build and validation

The requested profile is explicit in the build configuration:

```sh
cmake -S two-round-online -B two-round-online/build/fkos \
  -Demp-ot_DIR="$PWD/two-round-online/build/deps/emp-ot" \
  -DEMP_AG_OT_BACKEND=Ferret -DEMP_AG_TRIPLE_BACKEND=FKOS \
  -DEMP_AG_GMW_TRIPLE_BATCH_SIZE=350000 -DMLDSA_PREPROCESS_SSP=40
cmake --build two-round-online/build/fkos -j 4
ctest --test-dir two-round-online/build/fkos --output-on-failure
```

For a matched comparison, configure `HalfGate` with the same 40-bit
parameter and 350,000 maximum output batch. The historical default
80-bit/small-batch profile remains available and must not be presented
as a matched comparison against the 40-bit large-batch profile.

Tests cover every bank triple, small-circuit truth tables, two and seven
parties, random-bit/triple interleaving on one COT stream, invalid products,
invalid input bits, corrupted cut/sacrifice/combine openings, and use
after finish/abort. The full network test includes accepting and rejecting
MLDSA44 traces through actual preprocessing, saved WRK files, and the
unchanged two-flight online protocol. Public context agreement binds the
OT backend, triple backend, statistical parameter, and maximum batch size.
An additional exact-rational test checks 8,034 small-population correctness
and leakage probabilities against the parameter bounds, including the
production size-3 configuration and the insufficient 250,000-output case.

## Measured seven-party MLDSA44 comparison

Local arm64 Release builds, 2026-09-07; both use Ferret, a 40-bit
per-batch statistical parameter, and a 350,000 output-batch limit.
Communication below is aggregate application bytes **sent once across
all seven parties**, not sent plus received; aggregate received bytes
match exactly. Free edaBit provisioning and keygen are excluded.

| Provider | Before `c` bytes | After `c` bytes | Total bytes |
|---|---:|---:|---:|
| HalfGate | 511,308,129 | 3,162,828 | 514,470,957 |
| CutChoose, size 3+3 | 241,787,661 | 3,162,828 | 244,950,489 |

Before-`c` communication falls by **269,520,468 bytes (52.71%)**.
Both accepting and rejecting traces produce these counts, with exactly
1,024 public `w1` coefficients and 18,433 online output bits checked
against the plaintext specification. Online evaluation still has two
message flights, 54,271 garbled ANDs, and the same masked-output policy.
The Boolean work remains 251,904 Phase-1 ANDs (AND depth 80), 16,384
key-BtoA ANDs, and 54,271 WRK mask-product ANDs.

CutChoose before-`c` decomposition (disjoint rows):

| Work | Bytes sent |
|---|---:|
| COT setup | 10,396,281 |
| Random triple generation, including its checks | 75,627,552 |
| Phase-1 GMW evaluation | 2,863,392 |
| Key BtoA | 293,370 |
| WRK authenticated masks | 803,040 |
| WRK mask-product GMW evaluation | 572,586 |
| WRK offline delivery | 150,721,896 |
| Field opening, `w1`, input binding, context/finish/ready checks | 509,544 |
| **Total before `c`** | **241,787,661** |

The 75,627,552-byte triple-generation subtotal consists of 11,015,088
bytes for candidate COT authentication/checks, 15,241,254 for hash-product
corrections, 15,480,486 for product authentication, 10,080 for permutation
coins, 3,696 for the random cut, 30,487,338 for correctness sacrifice,
and 3,389,610 for leakage combination. It uses 8,754,841 explicitly drawn
authentication COTs per directed pair, including 45,568 discarded check
padding bits. These figures are subdivisions, not additional costs.

HalfGate triple generation sends 345,156,756 bytes, so the triple-stage
saving is **78.09%**. In the accepting local run, the maximum party's
triple-stage timer is 46.45 seconds for CutChoose versus 7.80 seconds for
HalfGate. The maximum per-party sum of all timed preprocessing stages
is 50.92 versus 12.02 seconds. These are single-run localhost timings,
not WAN benchmarks; the new provider reduces traffic but increases local
computation and peak memory. Party 1 (the evaluator) sends 12,989,538
before-`c` bytes and receives 163,751,148; the other parties each send
38,116,473–38,149,568. Aggregate traffic must not be mistaken for traffic
per party.

The final CutChoose build passed all nine non-preprocessing CTests and
the complete preprocessing network test separately. The matched HalfGate
network test also passed accepting/rejecting traces and both injected
authentication faults; mixed triple backends abort at context agreement
before COT setup, without persisting online state. Local raw measurements
remain under `two-round-online/build/fkos/preprocessing-network-0scpeajb`
and `two-round-online/build/halfgate40/preprocessing-network-qi68145q`.
