set(gfalib_srcs gfakluge.cpp)
add_library(gfakluge STATIC ${gfalib_srcs})

add_executable(pufferfish-index PufferFishIndex.cpp FastxParser.cpp)
target_link_libraries(pufferfish-index z sdsl gfakluge)

add_executable(test-gfa OurGFAReader.cpp)

