#!/bin/sh -e

myprefix=$HOME/local
PFFT_VERSION=1.0.6-alpha
FFTW_VERSION=3.3.3
INSTDIR=$myprefix/pfft-$PFFT_VERSION
FFTWDIR=$myprefix/fftw-$FFTW_VERSION
TMP="tmp-pfft-$PFFT_VERSION"

# bash check if directory exists
if [ -d $TMP ]; then
        echo "Directory $TMP exists. Cleanup before new build."
        exit 1
fi

mkdir $TMP && cd $TMP
cd ../.. && ./bootstrap.sh && cd -
../../configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR --disable-shared FC=mpif90 CC=mpicc MPICC=mpicc MPIFC=mpif90

make -j 4
make install
# make check
# ./configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR FC="mpif90 -f90=gfortran" CC="mpicc -cc=gcc"&& make && make check && make install
