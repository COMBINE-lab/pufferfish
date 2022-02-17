project(twopaco CXX)
cmake_minimum_required(VERSION 3.9)

set(CMAKE_PROJECT_NAME twopaco)

if(${TBB_FOUND})
  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  message("in pufferfish::twopaco::graphconstructor, TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
endif()

include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIRS} "../common")

add_library(twopaco STATIC ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
install(TARGETS twopaco DESTINATION lib EXPORT twopaco-targets)
install(EXPORT twopaco-targets DESTINATION lib/twopaco)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")
