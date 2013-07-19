#!/bin/sh -e

myprefix=$HOME/local
PFFT_VERSION=1.0.6-alpha
FFTW_VERSION=3.3.3
INSTDIR=$myprefix/pfft-$PFFT_VERSION-debug
FFTWDIR=$myprefix/fftw-$FFTW_VERSION-debug
TMP="tmp-pfft-$PFFT_VERSION-debug"

# bash check if directory exists
if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
	read answer
	if [ ${answer} = "y" ]; then
		rm -rf $TMP
	else
		echo "Program aborted."
		exit 1
	fi
fi

mkdir $TMP && cd $TMP
cd ../../ && ./bootstrap.sh && cd -

../../configure --build=powerpc64-bgp-linux-gnu --host=powerpc-ibm-cnk \
--prefix=$INSTDIR --with-fftw3=$FFTWDIR --enable-debug --disable-shared \
MPICC=mpixlc_r MPIFC=mpixlf90_r MPICXX=mpixlcxx_r \
CFLAGS='-O0 -g -qfullpath -qdbxextra -qcheck' FCFLAGS='-O0 -g -qfullpath'

make -j 4
make install

