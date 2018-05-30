# image: COMBINE-lab/salmon
# This dockerfile is based on the one created by
# Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile)
#FROM insilicodb/bioconductor 
FROM ubuntu:16.04

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf
ENV SALMON_VERSION 0.9.0

# pufferfish executable will be available in /home/pufferfish/build/src
# salmon binary will be installed in /home/salmon/bin/salmon

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

# SeqLib, the c++ wrapper for htslib
RUN git clone --recursive https://github.com/walaj/SeqLib.git && \
	cd SeqLib && \
	./configure && \
	make CXXFLAGS='-std=c++11' && \ 
	make install && \
	make seqtools 

# Pufferfish
RUN git clone https://github.com/COMBINE-lab/pufferfish.git && \
	cd pufferfish && \
	git checkout skoufos && \
	git clone https://github.com/simongog/sdsl-lite.git && \
	cd sdsl-lite && \ 
	./install.sh ../ && \
	cd .. && \
	mkdir build && \
	cd build && \
	cmake ../ && \
	make

# Salmon
RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

ENV PATH home/salmon-${SALMON_VERSION}/bin:${PATH}
