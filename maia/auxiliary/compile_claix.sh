#!/bin/bash

cd ~

ROOT_DIR=${PWD}

if [ "${1}" = "gnu" ]; then
    echo "Using GNU compiler"
    module unload intel intelmpi gcc openmpi
    module load foss/2022a
    compiler="gnu"
else
    echo "Using Intel compiler"
    module unload intel intelmpi gcc openmpi
    module load intel/2022a
    module load GCCcore/.11.3.0
    compiler="intel"
fi


if [ ! -d "${ROOT_DIR}/libraries/"]
then
    mkdir "${ROOT_DIR}/libraries/"
    echo "Created ${ROOT_DIR}/libraries/"
else 
    echo "${ROOT_DIR}/libraries/ exists"
fi

export FFTW_INST="${ROOT_DIR}/libraries/fftw_${compiler}"
echo "Installing FFTW in ${FFTW_INST}"

echo ${FFTW_INST}

export CXX=g++ 

wget http://www.fftw.org/fftw-3.3.7.tar.gz

tar xzf fftw-3.3.7.tar.gz

cd fftw-3.3.7

./configure --prefix="${FFTW_INST}" --enable-mpi --with-mpi=mpigxx

make -j

make install

cd ..

 
