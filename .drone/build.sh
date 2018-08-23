#!/bin/bash
source /hbb_exe/activate

set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"

# SeqLib, the c++ wrapper for htslib
echo "[Drone build] maping seqLib"
mkdir -p external_dependencies
cd external_dependencies
git clone --recursive https://github.com/walaj/SeqLib.git
cd SeqLib
./configure
make CXXFLAGS='-std=c++11'
make install
make seqtools
cd ../..

echo "[Drone build] Installikng sdsl"
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh ../
cd ..

echo "[Drone build] making build directory"
mkdir build
cd build

echo "[Drone build] cmake configuration"
cmake ..

echo "[Drone build] making pufferfish and installing locally (this could take a while)"
make
