#!/bin/bash
source /hbb_exe/activate

set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"

echo "[Drone build] Installikng sdsl"
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh ../
cd ..

echo "[Drone build] making build directory"
mkdir build
cd build

echo "[Drone build] cmake configuration"
cmake -DDO_QUIET_MAKE=TRUE ..

echo "[Drone build] making pufferfish and installing locally (this could take a while)"
make install
