cmake_minimum_required(VERSION 3.9)
#cmake_policy(SET CMP0069 NEW)

if(DEFINED ENV{CC})
  set(CC $ENV{CC})
else()
  set(CC gcc)
endif()
message("CC: ${CC}")

set(CC_VERSION "")
if(${CC} MATCHES ^gcc-)
  string(REGEX REPLACE "gcc-" "" CC_VERSION ${CC})
endif()
message("CC version: ${CC_VERSION}")


project(pufferfish)
include(ExternalProject)

# auto-populate version:
# from https://stackoverflow.com/questions/47066115/cmake-get-version-from-multi-line-text-file
file(READ "current_version.txt" ver)

string(REGEX MATCH "VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})

set(CPACK_PACKAGE_VERSION_MAJOR ${ver_major})
set(CPACK_PACKAGE_VERSION_MINOR ${ver_minor})
set(CPACK_PACKAGE_VERSION_PATCH ${ver_patch})

set(CPACK_PACKAGE_VERSION "${ver_major}.${ver_minor}.${ver_patch}")
message("version: ${CPACK_PACKAGE_VERSION}")


set(PROJECT_VERSION ${CPACK_PACKAGE_VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Maryland")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pufferfish — a tool for indexing and alignment based on the compacted colored de Bruijn Graph")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE "Release")
endif()

# We require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc -static-libstdc++")
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if (BUILD_PUFF_FOR_SALMON) 
    message("Building basic pufferfish components for salmon")
endif()

if (NOT GAT_SOURCE_DIR)
    message("Top-level source directory variable not set externally; setting it to ${CMAKE_CURRENT_SOURCE_DIR}")
    set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set(KSW_FLAGS "-DHAVE_KALLOC")
#set(PF_CPP_FLAGS "-msse4;-W;-Wall;-Wextra;-Wpointer-arith;-Wunused;-Wwrite-strings;-Wno-unknown-pragmas;-Wno-unused-function")
set(PF_CPP_FLAGS "-W;-Wall;-Wextra;-Wpointer-arith;-Wunused;-Wwrite-strings;-Wno-unknown-pragmas;-Wno-unused-function")

if (BUILD_PUFF_FOR_SALMON)
    #list(APPEND PF_CPP_FLAGS "-Wno-unused-parameter") 
    #list(APPEND PF_CPP_FLAGS "-Wno-unused-variable") 
    list(APPEND PF_CPP_FLAGS "-Wno-reorder") 
    #list(APPEND PF_CPP_FLAGS "-Wno-sign-compare") 
    #list(APPEND PF_CPP_FLAGS "-Wno-sequence-point")
   list(APPEND PF_CPP_FLAGS "-DPUFFERFISH_SALMON_SUPPORT")
endif()

###
# Sanitizers BEGIN
###
if (ASAN_BUILD)
  list(APPEND PF_CPP_FLAGS "-fsanitize=address")
#list(APPEND PF_CPP_FLAGS "-fsanitize=undefined")
#set(CMAKE_LINK_FLAGS "-fsanitize=address")
#list(APPEND CMAKE_LINK_FLAGS "-fsanitize=undefined")
# list(APPEND PF_CPP_FLAGS "-g")
  set(ASAN_LIB "asan")
else()
  set(ASAN_LIB "")
endif()
###
# Sanitizers END
###

set(WARN_ALL_THINGS "-fdiagnostics-color=always;-Wall;-Wcast-align;-Wcast-qual;-Wconversion;-Wctor-dtor-privacy;-Wdisabled-optimization;-Wdouble-promotion;-Wextra;-Wformat=2;-Winit-self;-Wlogical-op;-Wmissing-declarations;-Wmissing-include-dirs;-Wno-sign-conversion;-Wnoexcept;-Wold-style-cast;-Woverloaded-virtual;-Wpedantic;-Wredundant-decls;-Wshadow;-Wstrict-aliasing=1;-Wstrict-null-sentinel;-Wstrict-overflow=5;-Wswitch-default;-Wundef;-Wno-unknown-pragmas;-Wuseless-cast;-Wno-unused-parameter")

#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wrestrict -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 

set(OPT_FLAGS "-DSTX_NO_STD_STRING_VIEW;-O3;-fPIC;-DNDEBUG;-funroll-loops;-ftree-vectorize;-fno-strict-aliasing")
set(DEBUG_FLAGS "-pg;-g;-gstabs")

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if(APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    list(APPEND TGT_COMPILE_FLAGS "-undefined dynamic_lookup;-Wno-unused-command-line-argument")
    set(LIBSALMON_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")
else()
    # We're in sane linux world
    set(SHARED_LIB_EXTENSION "so")
    set(LIBSALMON_LINKER_FLAGS "")
#    if (NOT BUILD_PUFF_FOR_SALMON)
#        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#        list(APPEND PF_CPP_FLAGS -static-libstdc++)
#        message("PF_CPP_FLAGS = ${PF_CPP_FLAGS}")
#    endif()
endif()

##
# Let us check the sha sum of our packages if we have the right tools
##
set(SHASUM ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_shasum.sh)


find_package(Jemalloc)
    if(JEMALLOC_FOUND)
      include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
endif()



if (NOT BUILD_PUFF_FOR_SALMON)
## Try and find TBB first
find_package(TBB 2018.0 COMPONENTS tbb tbbmalloc tbbmalloc_proxy)

##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if((NOT TBB_FOUND) OR (TBB_FOUND AND (TBB_VERSION VERSION_LESS 2018.0)))

    set(TBB_WILL_RECONFIGURE TRUE)
    # Set the appropriate compiler
    if(CLANG)
        set(TBB_COMPILER "clang")
    else()
        set(TBB_COMPILER "gcc")
    endif()

    message("Build system will fetch and build Intel Threading Building Blocks")
    message("==================================================================")
    # These are useful for the custom install step we'll do later
    set(TBB_SOURCE_DIR ${GAT_SOURCE_DIR}/external/oneTBB-2020.1)
    set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)

    if("${TBB_COMPILER}" STREQUAL "gcc")
        ## Don't know why it's a problem yet, but if we're using
        ## GCC, get rid of the DO_ITT_NOTIFY flag
        set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
    endif()

    set(TBB_CXXFLAGS "${TBB_CXXFLAGS} ${CXXSTDFLAG}")

    ExternalProject_Add(libtbb
            DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
            DOWNLOAD_COMMAND curl -k -L https://github.com/oneapi-src/oneTBB/archive/v2020.1.tar.gz -o v2020.1.tar.gz &&
            ${SHASUM} 7c96a150ed22bc3c6628bc3fef9ed475c00887b26d37bca61518d76a56510971 v2020.1.tar.gz &&
            tar -xzvf v2020.1.tar.gz
            ##
            #URL https://github.com/01org/tbb/archive/2018_U3.tar.gz
            #DOWNLOAD_NAME 2018_U3.tar.gz
            #URL_HASH SHA1=d6cf16a42ece60aad6a722b369e1a2aa753347b4
            #TLS_VERIFY FALSE
            ##
            SOURCE_DIR ${TBB_SOURCE_DIR}
            INSTALL_DIR ${TBB_INSTALL_DIR}
            PATCH_COMMAND "${TBB_PATCH_STEP}"
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make ${QUIET_MAKE} CXXFLAGS=${TBB_CXXFLAGS} lambdas=1 compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
            INSTALL_COMMAND sh -c "mkdir -p ${TBB_INSTALL_DIR}/include && mkdir -p ${TBB_INSTALL_DIR}/lib && cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
            BUILD_IN_SOURCE 1
            )

    set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
    ExternalProject_Add_Step(libtbb reconfigure
            COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
            DEPENDEES install
            )

    set(FETCHED_TBB TRUE)
endif()

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
    set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
    set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
    set(TBB_LIB_DIR ${TBB_INSTALL_DIR}/lib)
    #set(TBB_LIBRARIES tbb tbbmalloc)
    set(TBB_LIBRARIES ${TBB_INSTALL_DIR}/lib/libtbb.${SHARED_LIB_EXTENSION}
			${TBB_INSTALL_DIR}/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
	   )
    message("TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
    message("TBB_LIBRARY_DIRS = ${TBB_LIBRARY_DIRS}")
endif()

##
#  The libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
if(TBB_RECONFIGURE)
    unset(TBB_FOUND CACHE)
    unset(TBB_INSTALL_DIR CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${GAT_SOURCE_DIR}/external/install)
    set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)
    message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
    find_package(TBB 2018.0 COMPONENTS tbb tbbmalloc tbbmalloc_proxy)
endif()
message("TBB_LIBRARIES = ${TBB_LIBRARIES}")


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
BUILD_COMMAND     make CXXFLAGS='-std=c++14'
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
endif()

set(LIB_FLAGS "-L ${CMAKE_CURRENT_SOURCE_DIR}/external/lib")
set(PUFF_RELEASE_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${PF_CPP_FLAGS};${OPT_FLAGS}")
set(PUFF_DEBUG_FLAGS "${DEFINE_FLAGS};${KSW_FLAGS};${PF_CPP_FLAGS};${DEBUG_FLAGS}")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/ntcard)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/twopaco)

if (NOT BUILD_PUFF_FOR_SALMON) 
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/asyncplusplus)
endif()

include(CPack)

include_directories(include)
link_directories(lib)
add_subdirectory(src)
