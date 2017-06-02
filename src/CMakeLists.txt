set(gfalib_srcs gfakluge.cpp)
set(pufferfish_lib_srcs GFAPufferFishIndex.cpp FastxParser.cpp OurGFAReader.cpp PufferFS.cpp)
add_library(gfakluge STATIC ${gfalib_srcs})
add_library(puffer STATIC ${pufferfish_lib_srcs})

add_executable(pufferfish-index-old PufferFishIndex.cpp FastxParser.cpp)
add_executable(pufferfish Pufferfish.cpp)
add_executable(fixFasta FixFasta.cpp FastxParser.cpp xxhash.c)

target_link_libraries(pufferfish-index-old z sdsl)
target_link_libraries(pufferfish z puffer sdsl)
target_link_libraries(fixFasta z)

#add_executable(test-gfa OurGFAReader.cpp)

