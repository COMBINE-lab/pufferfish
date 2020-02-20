project(twopaco CXX)
cmake_minimum_required(VERSION 3.9)

set(CMAKE_PROJECT_NAME twopaco)

include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIRS} "../common")
link_directories(${TBB_LIB_DIR})
add_library(twopaco STATIC ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
#add_executable(twopaco ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
#target_link_libraries(twopaco  "tbb")
install(TARGETS twopaco DESTINATION lib EXPORT twopaco-targets)
install(EXPORT twopaco-targets DESTINATION lib/twopaco)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")
