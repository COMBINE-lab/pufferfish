add_executable(pufferfish-index PufferFishIndex.cpp FastxParser.cpp)
target_link_libraries(pufferfish-index z sdsl)

