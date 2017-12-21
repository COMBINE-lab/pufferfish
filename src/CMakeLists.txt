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
    PufferfishIndex.cpp 
    PufferfishSparseIndex.cpp 
    Util.cpp 
    PufferfishValidate.cpp 
    PufferfishTestLookup.cpp 
    FastxParser.cpp 
    OurGFAReader.cpp 
    PufferFS.cpp 
    xxhash.c 
    GFAConverter.cpp
    PufferfishAligner.cpp
    edlib.cpp
	RefSeqConstructor.cpp
	)

add_library(puffer STATIC ${pufferfish_lib_srcs})

#add_executable(pufferfish-index-old PufferFishIndexer.cpp FastxParser.cpp)
add_executable(bcalm_pufferize BCALMPufferizer.cpp FastxParser.cpp)
add_executable(kswcli cli.cpp)
add_executable(pufferfish Pufferfish.cpp)
add_executable(fixFasta FixFasta.cpp FastxParser.cpp xxhash.c)
add_executable(myGFAtester MyGFATester.cpp FastxParser.cpp)
#add_executable(myGraphtester MyGraphTester.cpp)
add_executable(pufferize Pufferizer.cpp)
add_executable(edgedensity EdgeDensity.cpp)

#target_link_libraries(pufferfish-index-old z sdsl)
target_link_libraries(pufferfish puffer sdsl ksw2 ${JEMALLOC_LIBRARIES} z)
target_link_libraries(bcalm_pufferize z)
target_link_libraries(fixFasta z)
target_link_libraries(myGFAtester z sdsl puffer)
target_link_libraries(pufferize puffer z)
target_link_libraries(edgedensity z puffer)
target_link_libraries(kswcli z ksw2)
#target_link_libraries(myGraphtester puffer)


#add_executable(test-gfa OurGFAReader.cpp)

