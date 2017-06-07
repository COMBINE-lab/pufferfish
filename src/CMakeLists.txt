#set(gfalib_srcs gfakluge.cpp)
#add_library(gfakluge STATIC ${gfalib_srcs})
set(pufferfish_lib_srcs GFAPufferFishIndexer.cpp PufferfishIndex.cpp Util.cpp PufferfishValidate.cpp FastxParser.cpp OurGFAReader.cpp PufferFS.cpp xxhash.c)
add_library(puffer STATIC ${pufferfish_lib_srcs})

#add_executable(pufferfish-index-old PufferFishIndexer.cpp FastxParser.cpp)
add_executable(pufferfish Pufferfish.cpp)
add_executable(fixFasta FixFasta.cpp FastxParser.cpp xxhash.c)
add_executable(myGFAtester MyGFATester.cpp FastxParser.cpp)
add_executable(myGraphtester MyGraphTester.cpp)
add_executable(pufferizer k2k-1Converter.cpp)

#target_link_libraries(pufferfish-index-old z sdsl)
target_link_libraries(pufferfish z puffer sdsl)
target_link_libraries(fixFasta z)
target_link_libraries(myGFAtester z sdsl puffer)
target_link_libraries(pufferizer z puffer)
target_link_libraries(myGraphtester puffer)


#add_executable(test-gfa OurGFAReader.cpp)

