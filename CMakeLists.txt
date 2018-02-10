cmake_minimum_required(VERSION 3.1)

project(pufferfish)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# We require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(KSW_FLAGS "-march=native -DHAVE_KALLOC")
set(RBF_CPP_FLAGS "-g -mbmi2 -msse4 -pthread -std=c++11 -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -openmp -Wno-unknown-pragmas -Wno-unused-function")
set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 

#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wrestrict -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 


set(OPT_FLAGS "-O3 -ffast-math -fPIC -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
#set(OPT_FLAGS "-fprofile-generate -O3 -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
#set(OPT_FLAGS "-fprofile-use -fprofile-correction -O3 -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
set(DEBUG_FLAGS "-pg -g -gstabs")

set(RBF_CPP_FLAGS "${KSW_FLAGS} ${RBF_CPP_FLAGS} ${OPT_FLAGS}")

set(CMAKE_CXX_FLAGS ${RBF_CPP_FLAGS})

find_package(JeMalloc)
    if(JEMALLOC_FOUND)
      include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
endif()

include_directories(include)
link_directories(lib)
add_subdirectory(src)
