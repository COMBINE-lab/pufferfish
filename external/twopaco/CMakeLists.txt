cmake_minimum_required(VERSION 3.9)
project(root)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE RELEASE)
add_subdirectory(graphdump)
add_subdirectory(graphconstructor)
include(CPack)

#install(
#  FILES
#  ${CMAKE_BINARY_DIR}/graphconstructor/libtwopaco.a
#  ${CMAKE_BINARY_DIR}/graphdump/libgraphdump.a
#  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
