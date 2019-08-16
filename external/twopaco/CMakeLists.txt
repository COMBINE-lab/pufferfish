cmake_minimum_required(VERSION 3.9)
project(root)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE RELEASE)
add_subdirectory(graphdump)
add_subdirectory(graphconstructor)
include(CPack)

if (BUILD_PUFF_FOR_SALMON)
add_dependencies(graphdump libtbb)
add_dependencies(twopaco libtbb)
endif()

#install(
#  FILES
#  ${CMAKE_BINARY_DIR}/graphconstructor/libtwopaco.a
#  ${CMAKE_BINARY_DIR}/graphdump/libgraphdump.a
#  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
