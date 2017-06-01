set(gfalib_srcs gfakluge.cpp)
add_library(gfakluge STATIC ${gfalib_srcs})

add_executable(pufferfish-index-old PufferFishIndex.cpp FastxParser.cpp)
add_executable(pufferfish-index GFAPufferFishIndex.cpp FastxParser.cpp OurGFAReader.cpp)
target_link_libraries(pufferfish-index-old z sdsl)
target_link_libraries(pufferfish-index z sdsl)

#add_executable(test-gfa OurGFAReader.cpp)

