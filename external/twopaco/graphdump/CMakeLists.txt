project(graphdump CXX)
cmake_minimum_required(VERSION 3.9)

set(CMAKE_PROJECT_NAME graphdump)

include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIRS} "../common" "./" )
add_library(graphdump STATIC graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)

install(TARGETS graphdump DESTINATION lib EXPORT graphdump-targets)
install(EXPORT graphdump-targets DESTINATION lib/graphdump)
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")

if (FETCHED_TBB)
  add_dependencies(graphdump libtbb)
else()
  add_dependencies(graphdump TBB::tbb)
  add_dependencies(graphdump TBB::tbbmalloc)
endif()
