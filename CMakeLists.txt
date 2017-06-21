cmake_minimum_required(VERSION 3.0)

project(pufferfish)

set(RBF_CPP_FLAGS "-pthread -std=c++11 -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -openmp -Wno-unknown-pragmas")
set(OPT_FLAGS "-O3 -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
set(DEBUG_FLAGS "-pg -g -gstabs")

set(RBF_CPP_FLAGS "${RBF_CPP_FLAGS} ${OPT_FLAGS}")

set(CMAKE_CXX_FLAGS ${RBF_CPP_FLAGS})

include_directories(include)
link_directories(lib)
add_subdirectory(src)
