/* Compile the xxHash implementation from SSHash's newer version.
 * XXH_STATIC_LINKING_ONLY exposes internal structs needed by XXH_IMPLEMENTATION. */
#define XXH_STATIC_LINKING_ONLY
#define XXH_IMPLEMENTATION
#include "xxhash.h"
