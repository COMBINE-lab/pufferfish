#!/bin/bash
source /hbb_shlib/activate

set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"


echo "[Drone build] Installing ThreePaCo"
yum install tbb
git clone https://github.com/fataltes/TwoPaCo.git
cd TwoPaCo
mkdir build
cd build
cmake ../src
make
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
