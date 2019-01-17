cmake_minimum_required(VERSION 3.9)
#cmake_policy(SET CMP0069 NEW)

project(pufferfish)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE "Release")
endif()

# We require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(KSW_FLAGS "-march=native;-DHAVE_KALLOC")
set(RBF_CPP_FLAGS "-g;-mbmi2;-msse4;-W;-Wall;-Wextra;-Wpointer-arith;-Wunused;-Wwrite-strings;-openmp;-Wno-unknown-pragmas;-Wno-unused-function")
set(WARN_ALL_THINGS "-fdiagnostics-color=always;-Wall;-Wcast-align;-Wcast-qual;-Wconversion;-Wctor-dtor-privacy;-Wdisabled-optimization;-Wdouble-promotion;-Wextra;-Wformat=2;-Winit-self;-Wlogical-op;-Wmissing-declarations;-Wmissing-include-dirs;-Wno-sign-conversion;-Wnoexcept;-Wold-style-cast;-Woverloaded-virtual;-Wpedantic;-Wredundant-decls;-Wshadow;-Wstrict-aliasing=1;-Wstrict-null-sentinel;-Wstrict-overflow=5;-Wswitch-default;-Wundef;-Wno-unknown-pragmas;-Wuseless-cast;-Wno-unused-parameter")

#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wrestrict -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 

set(OPT_FLAGS "-O3;-ffast-math;-fPIC;-DNDEBUG;-funroll-loops;-mmmx;-msse;-msse2;-msse3;-msse4;-msse4.2;-march=native;-fno-strict-aliasing")
set(DEBUG_FLAGS "-pg;-g;-gstabs")

find_package(Jemalloc)
    if(JEMALLOC_FOUND)
      include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
endif()

include(ExternalProject)
ExternalProject_Add(libseqlib
GIT_REPOSITORY https://github.com/COMBINE-lab/SeqLib.git
GIT_TAG        master
UPDATE_COMMAND ""
UPDATE_DISCONNECTED 1
BUILD_IN_SOURCE TRUE
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/seqlib
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/seqlib
INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
CONFIGURE_COMMAND ./configure 
BUILD_COMMAND     make -j8 CXXFLAGS='-std=c++14'
INSTALL_COMMAND   mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include && cp src/libseqlib.a <INSTALL_DIR>/lib && 
                  cp htslib/libhts.a <INSTALL_DIR>/lib &&
                  cp -r SeqLib <INSTALL_DIR>/include &&
                  cp -r json <INSTALL_DIR>/include &&
                  cp -r htslib <INSTALL_DIR>/include
)
ExternalProject_Add(libSetCover
        GIT_REPOSITORY https://github.com/martin-steinegger/setcover.git
        GIT_TAG        master
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        UPDATE_DISCONNECTED 1
        BUILD_IN_SOURCE TRUE
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/setcover
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/setcover
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        BUILD_COMMAND     make -j8
        INSTALL_COMMAND     ar rcs libSetCover.a SetCover.o &&
        mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include &&
        cp libSetCover.a <INSTALL_DIR>/lib &&
        cp SetCover.h <INSTALL_DIR>/include &&
        cp SetElement.h <INSTALL_DIR>/include &&
        cp LinearMultiArray.h <INSTALL_DIR>/include
        )
set(LIB_FLAGS "-L ${CMAKE_CURRENT_SOURCE_DIR}/external/lib")
set(PUFF_RELEASE_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${RBF_CPP_FLAGS};${OPT_FLAGS}")
set(PUFF_DEBUG_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${RBF_CPP_FLAGS};${DEBUG_FLAGS}")

include_directories(include)
link_directories(lib)
add_subdirectory(src)
