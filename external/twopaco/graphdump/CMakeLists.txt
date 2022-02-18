project(graphdump CXX)
cmake_minimum_required(VERSION 3.9)

set(CMAKE_PROJECT_NAME graphdump)

#if(${TBB_TARGET_EXISTED})
#  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
#  message("in pufferfish::twopaco::graphdump, TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
#endif()

include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIRS} "../common" "./" )
#link_directories(${TBB_LIB_DIR})

add_library(graphdump STATIC graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)

#add_executable(graphdump graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)
#target_link_libraries(graphdump  "tbb")
install(TARGETS graphdump DESTINATION lib EXPORT graphdump-targets)
install(EXPORT graphdump-targets DESTINATION lib/graphdump)
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")

add_dependencies(graphdump libtbb-install)
