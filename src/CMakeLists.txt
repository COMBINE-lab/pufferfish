
####
# Deprecated
####
#INCLUDE_DIRECTORIES(/home/fatemeh/others_projects/SeqLib/ /home/fatemeh/others_projects/SeqLib/htslib/ )
#SET(SEQLIB_DIR /home/rob/nemo/external/SeqLib/bin)
#INCLUDE_DIRECTORIES(../external/include/SeqLib/)
#SET(SEQLIB_DIR ../external_dependencies/SeqLib/bin)
#set(gfalib_srcs gfakluge.cpp)
#add_library(gfakluge STATIC ${gfalib_srcs})
####
# End Deprecated
####

# check if we know how to do IPO
include(CheckIPOSupported)
check_ipo_supported(RESULT HAS_IPO)

link_directories(
		${CMAKE_CURRENT_SOURCE_DIR}/../external/install/lib
		${TBB_LIBRARIES}
)

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

set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_FLAGS "-O3 -msse2 -mno-sse4.1")
set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;KSW_SSE2_ONLY;HAVE_KALLOC")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_FLAGS "-O3 -msse4.1")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
set_target_properties(ksw2pp_basic PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")

set_target_properties(ksw2pp_basic PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include)
set_target_properties(ksw2pp_sse4 PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include) 

# Build the ksw2pp library
add_library(ksw2pp STATIC $<TARGET_OBJECTS:ksw2pp_sse2> $<TARGET_OBJECTS:ksw2pp_sse4> $<TARGET_OBJECTS:ksw2pp_basic>)
set_target_properties(ksw2pp PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
if(HAS_IPO)
  set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

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
    PufferfishTestLookup.cpp 
    FastxParser.cpp 
#    PufferfishGFAReader.cpp
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

add_library(puffer STATIC ${pufferfish_lib_srcs})


target_compile_options(puffer PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(puffer PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(bcalm_pufferize BCALMPufferizer.cpp FastxParser.cpp)
target_compile_options(bcalm_pufferize PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(bcalm_pufferize PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

#add_executable(kswcli cli.cpp)

add_executable(pufferfish Pufferfish.cpp)
target_compile_options(pufferfish PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(pufferfish PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

## Apparently, CMake has a multi-line comment format.  I didn't know this but
## CLion did (@fataltes).
#[[
add_executable(fixFasta FixFasta.cpp FastxParser.cpp xxhash.c)
target_compile_options(fixFasta PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(fixFasta PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
]]

add_executable(edgedensity EdgeDensity.cpp)
target_compile_options(edgedensity PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(edgedensity PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(edgedensity2 EdgeDensityCustom.cpp)
target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(krakmap KrakMap.cpp Taxa.cpp)
target_compile_options(krakmap PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(krakmap PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")


add_executable(cedar TargetGroup.cpp Cedar.cpp Taxa.cpp)
target_include_directories(cedar PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/install/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/install/include/htslib
	${TBB_INCLUDE_DIRS}
)
target_compile_options(cedar PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(cedar PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(filtersam testing/SAMFilter.cpp)
target_include_directories(filtersam PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/install/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/install/include/htslib
    )
target_compile_options(filtersam PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(filtersam PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
add_dependencies(filtersam libseqlib)

#get_target_property(dirs filtersam INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()

#add_executable(readKrk readKrk.cpp taxa.cpp)

add_executable(rank_test rank_test.cpp rank9b.cpp rank9sel.cpp)
target_compile_options(rank_test PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(rank_test PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")
target_link_libraries(rank_test)

target_link_libraries(pufferfish Threads::Threads puffer ksw2pp ${JEMALLOC_LIBRARIES} z twopaco graphdump ntcard ${TBB_LIBRARIES})
target_link_libraries(bcalm_pufferize Threads::Threads z)
#target_link_libraries(fixFasta Threads::Threads z)
target_link_libraries(edgedensity Threads::Threads z puffer)
target_link_libraries(edgedensity2 Threads::Threads z puffer)
target_link_libraries(krakmap Threads::Threads puffer)
target_link_libraries(cedar ${TBB_LIBRARIES} Threads::Threads seqlib hts puffer bz2 lzma z SetCover)
target_link_libraries(filtersam Threads::Threads seqlib hts puffer bz2 lzma z )

