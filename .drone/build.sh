#!/bin/bash
set -e

source /hbb_exe/activate

set -x

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"
echo "[Drone build] making build directory"

yum install -y --quiet wget
wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
rpm -i --quiet epel-release-5-4.noarch.rpm

mkdir build
cd build

echo "[Drone build] cmake configuration"

cmake -DDO_QUIET_MAKE=TRUE ..

echo "[Drone build] making pufferfish and installing locally (this could take a while)"

