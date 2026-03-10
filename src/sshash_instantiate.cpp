// Explicit template instantiations for SSHash dictionary types.
// When SSHash is used as a library (not as a standalone tool),
// these instantiations are needed so that the linker can find
// the template method implementations defined in .cpp files.

#include "dictionary_types.hpp"

// Re-include the template method definitions so they're available
// for instantiation in this translation unit.
#include "dictionary.cpp"
#include "info.cpp"
#include "build.cpp"

namespace sshash {
template struct dictionary<default_kmer_t, decoded_offsets>;
}  // namespace sshash
