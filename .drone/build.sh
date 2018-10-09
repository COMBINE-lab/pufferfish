#!/bin/bash
set -e

source /hbb_exe/activate

set -x

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"


echo "[Drone build] Installing ThreePaCo"
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

wget https://github.com/01org/tbb/releases/download/2019_U1/tbb2019_20181003oss_lin.tgz
tar -xvzf tbb2019_20181003oss_lin.tgz
export LD_LIBRARY_PATH=${PWD}/tbb2019_20181003oss/lib/intel64/gcc4.7/
git clone https://github.com/fataltes/TwoPaCo.git
cd TwoPaCo
git checkout pufferize
mkdir build
cd build
cmake ../src/ -DTBB_LIB_DIR=/tbb2019_20181003oss/lib/intel64/gcc4.7/ -DTBB_INCLUDE_DIR=/tbb2019_20181003oss/include/
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
