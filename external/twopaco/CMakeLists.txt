cmake_minimum_required(VERSION 3.9)
project(root)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RELEASE)
add_subdirectory(graphdump)
add_subdirectory(graphconstructor)
include(CPack)

if (BUILD_PUFF_FOR_SALMON)
  if (${TBB_TARGET_EXISTED})
  add_dependencies(graphdump tbb)
  add_dependencies(twopaco tbb)
  else ()
  message("adding externally-fetched tbb as dependency of twopaco")
  add_dependencies(graphdump libtbb)
  add_dependencies(twopaco libtbb)
  endif()
else()
  if (${TBB_TARGET_EXISTED})
  add_dependencies(graphdump libtbb_rescan)
  add_dependencies(twopaco libtbb_rescan)
  else ()
  message("adding externally-fetched tbb as dependency of twopaco")
  add_dependencies(graphdump libtbb_rescan)
  add_dependencies(twopaco libtbb_rescan)
  endif()
endif()


#install(
#  FILES
#  ${CMAKE_BINARY_DIR}/graphconstructor/libtwopaco.a
#  ${CMAKE_BINARY_DIR}/graphdump/libgraphdump.a
#  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
