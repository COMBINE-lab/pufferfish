# check if we know how to do IPO
include(CheckIPOSupported)
check_ipo_supported(RESULT HAS_IPO)

function(pufferfish_configure_target target_name)
  target_link_libraries(${target_name}
    PRIVATE
      pufferfish_project_options
      pufferfish_repo_headers
  )
endfunction()

function(pufferfish_configure_standalone_target target_name)
  pufferfish_configure_target(${target_name})
  target_link_libraries(${target_name}
    PRIVATE
      pufferfish_external_link_dirs
  )
endfunction()

set (KSW2PP_BASIC_LIB_SRCS
ksw2pp/kalloc.c
ksw2pp/ksw2_extd.c
ksw2pp/ksw2_extz.c
ksw2pp/ksw2_gg.c
ksw2pp/ksw2_gg2.c
ksw2pp/ksw2_gg2_sse.c
ksw2pp/KSW2Aligner.cpp
)

set (KSW2PP_ADVANCED_LIB_SRCS
ksw2pp/ksw2_extd2_sse.c
ksw2pp/ksw2_extf2_sse.c
ksw2pp/ksw2_extz2_sse.c
)

add_library(ksw2pp_sse2 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_sse4 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_basic OBJECT ${KSW2PP_BASIC_LIB_SRCS})

if (USE_ARM)
    target_compile_options(ksw2pp_sse2 PRIVATE -O3)
    target_compile_definitions(ksw2pp_sse2 PRIVATE __SSE2__ KSW_CPU_DISPATCH KSW_USE_ARM KSW_SSE2_ONLY HAVE_KALLOC)
    target_compile_options(ksw2pp_sse4 PRIVATE -O3)
    target_compile_definitions(ksw2pp_sse4 PRIVATE __SSE2__ __SSE4_1__ KSW_CPU_DISPATCH KSW_USE_ARM HAVE_KALLOC)
    target_compile_definitions(ksw2pp_basic PRIVATE __SSE2__ KSW_CPU_DISPATCH KSW_USE_ARM HAVE_KALLOC)
else()
    target_compile_options(ksw2pp_sse2 PRIVATE -O3 -msse -msse2 -mno-sse4.1)
    target_compile_definitions(ksw2pp_sse2 PRIVATE KSW_CPU_DISPATCH KSW_SSE2_ONLY HAVE_KALLOC)
    target_compile_options(ksw2pp_sse4 PRIVATE -O3 -msse -msse2 -msse3 -mssse3 -msse4 -msse4.1)
    target_compile_definitions(ksw2pp_sse4 PRIVATE KSW_CPU_DISPATCH HAVE_KALLOC)
    target_compile_definitions(ksw2pp_basic PRIVATE KSW_CPU_DISPATCH HAVE_KALLOC)
endif()

foreach(ksw_target IN ITEMS ksw2pp_sse2 ksw2pp_sse4 ksw2pp_basic)
  pufferfish_configure_target(${ksw_target})
endforeach()
target_include_directories(ksw2pp_basic PRIVATE ${PUFFERFISH_ROOT_DIR}/include)
target_include_directories(ksw2pp_sse4 PRIVATE ${PUFFERFISH_ROOT_DIR}/include)
target_include_directories(ksw2pp_sse2 PRIVATE ${PUFFERFISH_ROOT_DIR}/include)

# Build the ksw2pp library
add_library(ksw2pp STATIC $<TARGET_OBJECTS:ksw2pp_sse2> $<TARGET_OBJECTS:ksw2pp_sse4> $<TARGET_OBJECTS:ksw2pp_basic>)
pufferfish_configure_target(ksw2pp)
target_compile_definitions(ksw2pp PRIVATE KSW_CPU_DISPATCH HAVE_KALLOC)
if(HAS_IPO)
  #set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()
target_include_directories(ksw2pp PUBLIC ${PUFFERFISH_ROOT_DIR}/include)
add_library(pufferfish::ksw2pp ALIAS ksw2pp)

set(pufferfish_lib_srcs 
    PufferfishIndexer.cpp 
    PufferfishBaseIndex.cpp
    PufferfishIndex.cpp 
    PufferfishSparseIndex.cpp 
    PufferfishLossyIndex.cpp
    edlib.cpp
    Util.cpp
		rank9sel.cpp
    rank9b.cpp
    PufferfishValidate.cpp
	  PufferfishStats.cpp
    PufferfishTestLookup.cpp 
    PufferfishExamine.cpp
    PufferfishKmerQuery.cpp
    FastxParser.cpp 
	  PufferfishBinaryGFAReader.cpp
    PufferFS.cpp
    xxhash.c
    FixFasta.cpp
    MemCollector.cpp
	  MemChainer.cpp
	  PuffAligner.cpp
	  PufferfishAligner.cpp
	  RefSeqConstructor.cpp
	  metro/metrohash64.cpp
)

message("add library puffer : ${pufferfish_lib_srcs}")
add_library(puffer STATIC ${pufferfish_lib_srcs})
pufferfish_configure_target(puffer)
target_include_directories(puffer
  PUBLIC
    ${PUFFERFISH_ROOT_DIR}/include
)
target_link_libraries(puffer
  PUBLIC
    ksw2pp
    twopaco
    graphdump
    ntcard
    Threads::Threads
    TBB::tbb
    TBB::tbbmalloc
)
add_library(pufferfish::puffer ALIAS puffer)

# fetch the directory where compact vector lives
get_property(COMPACT_VECTOR_INCLUDE_PATH TARGET graphdump PROPERTY COMPACT_VECTOR_DIR)
message("adding ${COMPACT_VECTOR_INCLUDE_PATH} to include path for libpuffer")
target_include_directories(puffer PUBLIC ${COMPACT_VECTOR_INCLUDE_PATH})

if (FETCHED_TBB)
  add_dependencies(puffer libtbb)
endif()

if (NOT PUFFERFISH_EMBEDDED)
message("Building extra components")

# The actual executable need not be built for salmon
# Also, the `-static-libgcc` flag is a problem under OSX/clang.
add_executable(pufferfish Pufferfish.cpp)
pufferfish_configure_standalone_target(pufferfish)

if (FETCHED_TBB)
add_dependencies(pufferfish libtbb)
else()
add_dependencies(pufferfish TBB::tbb)
add_dependencies(pufferfish TBB::tbbmalloc)
endif()



#if (PUFFERFISH_EMBEDDED)
#   target_include_directories(pufferfish PUBLIC "${PUFFERFISH_ROOT_DIR}/include")
#endif()


add_dependencies(pufferfish puffer)
add_dependencies(pufferfish ksw2pp) 
add_dependencies(pufferfish twopaco)
add_dependencies(pufferfish graphdump)
add_dependencies(pufferfish ntcard)

if(TBB_RECONFIGURE OR TBB_TARGET_EXISTED)
  target_link_libraries(pufferfish 
    ksw2pp
    puffer 
    z 
    twopaco 
    graphdump 
    ntcard 
    ${LIBRT} 
    ${ASAN_LIB} 
    ${JEMALLOC_LIBRARIES} 
    ksw2pp
    puffer 
    ${CMAKE_DL_LIBS} 
    TBB::tbb
    TBB::tbbmalloc
    Threads::Threads 
  )
endif()

add_executable(bcalm_pufferize BCALMPufferizer.cpp FastxParser.cpp)
pufferfish_configure_standalone_target(bcalm_pufferize)

#add_executable(edgedensity EdgeDensity.cpp)
#target_compile_options(edgedensity PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(edgedensity PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#add_dependencies(edgedensity puffer)
#

#add_executable(edgedensity2 EdgeDensityCustom.cpp)
#target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#add_dependencies(edgedensity2 puffer)

add_executable(filtersam testing/SAMFilter.cpp)
pufferfish_configure_standalone_target(filtersam)
target_include_directories(filtersam PRIVATE
    ${PUFFERFISH_ROOT_DIR}/external/install/include
    ${PUFFERFISH_ROOT_DIR}/external/install/include/htslib
    )
add_dependencies(filtersam libseqlib)
#add_dependencies(filtersam puffer)


#add_executable(krakmap KrakMap.cpp cedar/Taxa.cpp)
#target_compile_options(krakmap PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(krakmap PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(cedar cedar/TargetGroup.cpp cedar/Cedar.cpp cedar/Taxa.cpp)
pufferfish_configure_standalone_target(cedar)

target_include_directories(cedar PRIVATE
    ${PUFFERFISH_ROOT_DIR}/external/install/include
    ${PUFFERFISH_ROOT_DIR}/external/install/include/htslib
)
target_link_libraries(cedar Async++ Threads::Threads ${SEQLIB} ${SEQLIB_HTSLIB} puffer bz2 lzma z SetCover ${ASAN_LIB})
add_dependencies(cedar libseqlib)
add_dependencies(cedar libsetcover)
#add_dependencies(ceder Async++)
#add_dependencies(ceder seqlib)

#get_target_property(dirs filtersam INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()

#add_executable(readKrk readKrk.cpp taxa.cpp)


#[[
#add_executable(rank_test rank_test.cpp rank9b.cpp rank9sel.cpp)
#target_compile_options(rank_test PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
#target_compile_options(rank_test PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
#target_link_libraries(rank_test)
#]]

add_executable(getLineage constructLineage.cpp)
pufferfish_configure_standalone_target(getLineage)


target_link_libraries(bcalm_pufferize Threads::Threads z ${ASAN_LIB})
#target_link_libraries(fixFasta Threads::Threads z)


#if (${TBB_TARGET_EXISTED})
#target_link_libraries(edgedensity Threads::Threads z puffer ${ASAN_LIB})
#target_link_libraries(edgedensity2 Threads::Threads z puffer ${ASAN_LIB})
target_link_libraries(filtersam Threads::Threads ${SEQLIB} ${SEQLIB_HTSLIB} puffer bz2 lzma z ${ASAN_LIB})
add_dependencies(filtersam libseqlib)
#endif()


endif()
