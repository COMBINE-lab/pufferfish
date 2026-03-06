# Pufferfish For Salmon Integration Plan

## Goal

Expose a stable, modern CMake integration surface for Salmon without requiring:

- `scripts/fetchPufferfish.sh`
- copying headers and sources into `external/install`
- pufferfish-owned dependency bootstrap in the embedded build
- Salmon-specific hacks like injecting Salmon sources into pufferfish targets

## Current Embedded-Mode Problems

### Top-level CMake

File: `CMakeLists.txt`

- `BUILD_PUFF_FOR_SALMON` is a partial mode switch, not a clean embedded interface.
- `GAT_SOURCE_DIR` is used as a parent-project backchannel.
- `ExternalProject` bootstrap still owns TBB, liblzma, libbz2, SeqLib, and setcover in standalone mode.
- `find_package(Jemalloc)` is unconditional even though Salmon has moved to mimalloc.
- global `include_directories()` and `link_directories()` are still used.

### Source CMake

File: `src/CMakeLists.txt`

- `puffer` is the only library target and mixes indexing, alignment support, and Salmon-only concerns.
- `puffer` depends on the `graphdump` target property `COMPACT_VECTOR_DIR`.
- `BUILD_PUFF_FOR_SALMON` injects `${GAT_SOURCE_DIR}/src/LibraryFormat.cpp` into the pufferfish build.
- standalone executables and helper tools are gated in an ad hoc way rather than by target selection.

## Target End State

Provide these targets:

- `pufferfish::core`
- `pufferfish::align`
- `pufferfish::index`
- `pufferfish::ksw2pp`

With these rules:

- embedded builds do not fetch or install dependencies
- embedded builds consume imported targets provided by the parent
- include paths are target-scoped
- no parent-project source files are compiled into pufferfish targets
- standalone tools are optional and off by default for embedded use

## Proposed Build Options

- `PUFFERFISH_BUILD_TOOLS`
- `PUFFERFISH_BUILD_TESTS`
- `PUFFERFISH_BUILD_INDEXER`
- `PUFFERFISH_BUILD_ALIGN`
- `PUFFERFISH_EMBEDDED`
- `PUFFERFISH_USE_SYSTEM_DEPS`

Behavior:

- `PUFFERFISH_EMBEDDED=ON` disables dependency fetching, install rules, and standalone-only helpers.
- `PUFFERFISH_BUILD_TOOLS=OFF` in embedded mode.
- `PUFFERFISH_BUILD_TESTS=OFF` in embedded mode.
- `PUFFERFISH_BUILD_INDEXER=ON` and `PUFFERFISH_BUILD_ALIGN=ON` are left selectable by the parent.

## Mechanical Refactor Sequence

1. Replace global include/link state with interface targets.
2. Split `puffer` into smaller static libraries.
3. Remove `GAT_SOURCE_DIR` and `BUILD_PUFF_FOR_SALMON` source injection.
4. Replace `COMPACT_VECTOR_DIR` target property coupling with direct target include propagation.
5. Move standalone dependency bootstrap behind `if(NOT PUFFERFISH_EMBEDDED)`.
6. Export CMake package targets for standalone installation.

## Minimal Salmon-Facing Interface

Salmon currently needs these header families:

- `ProgOpts.hpp`
- `PufferfishBaseIndex.hpp`
- `PufferfishIndex.hpp`
- `PufferfishSparseIndex.hpp`
- `PufferfishLossyIndex.hpp`
- `Util.hpp`
- `SelectiveAlignmentUtils.hpp`
- `PuffAligner.hpp`
- `MemCollector.hpp`
- `MemChainer.hpp`
- `SAMWriter.hpp`
- `PufferfishConfig.hpp`
- `ksw2pp/KSW2Aligner.hpp`
- `metro/metrohash64.h`

Salmon currently needs these source families provided by pufferfish targets:

- index loading and query code
- MEM collection and chaining
- alignment formatting helpers
- `rank9b.cpp`
- `metrohash64.cpp`
- ksw2pp objects

Salmon should not need:

- copied source trees
- copied header trees
- standalone pufferfish executables
- SeqLib / cedar / filtersam / setcover

## First Code Changes To Make On This Branch

1. Introduce `PUFFERFISH_EMBEDDED` and alias it from `BUILD_PUFF_FOR_SALMON` temporarily.
2. Remove unconditional `find_package(Jemalloc)`.
3. Stop appending `${GAT_SOURCE_DIR}/src/LibraryFormat.cpp` in `src/CMakeLists.txt`.
4. Add target-scoped include directories for `ksw2pp` and `puffer`.
5. Restrict `ExternalProject` dependency setup to standalone mode only.
6. Export a clean `pufferfish::puffer` target first, then split it further.

## Salmon Migration Path

After the first pufferfish cleanup:

1. Salmon fetches a pinned pufferfish commit/tag via `FetchContent`.
2. Salmon links pufferfish targets directly.
3. Salmon removes `scripts/fetchPufferfish.sh`.
4. Salmon removes `external/install/include/pufferfish` and `external/install/src/pufferfish` from its dependency model.

