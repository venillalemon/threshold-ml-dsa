#ifndef EMP_AG_BACKEND_CRYPTO_FS_HASH_H__
#define EMP_AG_BACKEND_CRYPTO_FS_HASH_H__
//
// The hash behind every Fiat-Shamir / transcript-adjacent check digest in
// emp-ag: the F_eq window folds and their commit-open, the EchoBC view
// digests, the sampleRandom coin commitments, the shared-zero / d-vector /
// MAC-view digests. These are all "both parties compute and compare" hashes
// — an AGREEMENT-class choice with no interop pin beyond the peer — so they
// follow emp-tool's stack-wide EMP_FS_HASH selection (the same flag that
// picks the IOChannel transcript hash; PUBLIC define, inherited through
// emp-tool's exported target) instead of the general-purpose emp::Hash.
// With the flag unset this is exactly emp::Hash: byte-identical transcripts.
//
// General-purpose / non-check uses of emp::Hash (none in emp-ag today) and
// the MITCCRH garbling hash are deliberately NOT routed through this alias.
//
#include "emp-tool/runtime/crypto/hash.h"

// Building against an emp-tool that predates EMP_FS_HASH: follow the
// global default, which is what the old behavior was.
#ifndef EMP_FS_HASH_DEFAULT
#define EMP_FS_HASH_DEFAULT EMP_HASH_DEFAULT
#endif

namespace emp::ag {

using FsHash = emp::HashT<EMP_FS_HASH_DEFAULT>;

}  // namespace emp::ag
#endif  // EMP_AG_BACKEND_CRYPTO_FS_HASH_H__
