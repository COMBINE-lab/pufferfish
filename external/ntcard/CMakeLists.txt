cmake_minimum_required(VERSION 3.9)
project(ntcard)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE RELEASE)

include_directories("Common")

add_library(ntcard STATIC
        #Common/Uncompress.cpp
        #Common/Uncompress.h 
        Common/SignalHandler.cpp 
        Common/SignalHandler.h 
        Common/Fcontrol.cpp 
        Common/Fcontrol.h 
        Common/nthash.hpp 
        Common/ntHashIterator.hpp 
        nthll.cpp
        )

target_compile_options(ntcard PRIVATE -O3)

install(TARGETS ntcard DESTINATION lib EXPORT ntcard-targets)
install(EXPORT ntcard-targets DESTINATION lib/ntcard)

#if (NOT TBB_TARGET_EXISTED)
#  add_dependencies(ntcard libtbb)
#endif()