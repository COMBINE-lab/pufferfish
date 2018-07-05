
#INCLUDE_DIRECTORIES(/home/fatemeh/others_projects/SeqLib/ /home/fatemeh/others_projects/SeqLib/htslib/ )
LINK_DIRECTORIES(../external_dependencies/SeqLib/bin)
#SET(SEQLIB_DIR /home/rob/nemo/external/SeqLib/bin)
INCLUDE_DIRECTORIES(../external_dependencies/SeqLib/ ../external_dependencies/SeqLib/htslib/)
SET(SEQLIB_DIR ../external_dependencies/SeqLib/bin)
#set(gfalib_srcs gfakluge.cpp)
#add_library(gfakluge STATIC ${gfalib_srcs})
set(ksw2_lib_srcs
    kalloc.c
    KSW2Aligner.cpp
    ksw2_extd.c
    ksw2_extd2_sse.c
    ksw2_extf2_sse.c
    ksw2_exts2_sse.c
    ksw2_extz.c
    ksw2_extz2_sse.c
    ksw2_gg.c
    ksw2_gg2.c
    ksw2_gg2_sse.c
    )
add_library(ksw2 STATIC ${ksw2_lib_srcs})

set(pufferfish_lib_srcs 
    GFAPufferFishIndexer.cpp 
    PufferfishBaseIndex.cpp
    PufferfishIndex.cpp 
    PufferfishSparseIndex.cpp 
    PufferfishLossyIndex.cpp
    Util.cpp 
    PufferfishValidate.cpp 
    PufferfishTestLookup.cpp 
    FastxParser.cpp 
    OurGFAReader.cpp 
    PufferFS.cpp 
    xxhash.c 
    GFAConverter.cpp
    PufferfishAligner.cpp
    #edlib.cpp
	RefSeqConstructor.cpp
	Pufferize.cpp
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

add_executable(fixFasta FixFasta.cpp FastxParser.cpp xxhash.c)
target_compile_options(fixFasta PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(fixFasta PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(myGFAtester MyGFATester.cpp FastxParser.cpp)
target_compile_options(myGFAtester PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(myGFAtester PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

#add_executable(myGraphtester MyGraphTester.cpp)
add_executable(pufferize Pufferizer.cpp)
target_compile_options(pufferize PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(pufferize PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(edgedensity EdgeDensity.cpp)
target_compile_options(edgedensity PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(edgedensity PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(edgedensity2 EdgeDensityCustom.cpp)
target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(edgedensity2 PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

add_executable(krakmap krakMap.cpp taxa.cpp)
target_compile_options(krakmap PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(krakmap PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")


add_executable(cedar TargetGroup.cpp cedar.cpp taxa.cpp)
target_compile_options(cedar PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(cedar PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")


add_executable(filtersam testing/SAMFilter.cpp)
target_compile_options(filtersam PUBLIC "$<$<CONFIG:DEBUG>:${PUFF_DEBUG_FLAGS}>")
target_compile_options(filtersam PUBLIC "$<$<CONFIG:RELEASE>:${PUFF_RELEASE_FLAGS}>")

#add_executable(readKrk readKrk.cpp taxa.cpp)

#target_link_libraries(pufferfish-index-old z sdsl)
target_link_libraries(pufferfish Threads::Threads puffer sdsl ksw2 ${JEMALLOC_LIBRARIES} z)
target_link_libraries(bcalm_pufferize Threads::Threads z)
target_link_libraries(fixFasta Threads::Threads z)
target_link_libraries(myGFAtester Threads::Threads z sdsl puffer)
target_link_libraries(pufferize Threads::Threads puffer sdsl z)
target_link_libraries(edgedensity Threads::Threads z puffer)
target_link_libraries(edgedensity2 Threads::Threads z puffer)
#target_link_libraries(kswcli Threads::Threads z ksw2)
target_link_libraries(krakmap Threads::Threads puffer)
target_link_libraries(cedar Threads::Threads libseqlib.a libhts.a puffer bz2 lzma z )
target_link_libraries(filtersam Threads::Threads libseqlib.a libhts.a puffer bz2 lzma z )
#target_link_libraries(readKrk puffer bz2 lzma z ${SEQLIB_DIR}/libseqlib.a ${SEQLIB_DIR}/libhts.a)
#target_link_libraries(myGraphtester puffer)


#add_executable(test-gfa OurGFAReader.cpp)

