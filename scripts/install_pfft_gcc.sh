#!/bin/bash 

myprefix=$HOME/local
PFFT_VERSION=1.0.8-alpha
FFTW_VERSION=3.3.4

# ---- set path where the PFFT configure can be found ----------------- 
CONFDIR="$(pwd)/.."

# ---- set paths where the FFTW-MPI headers and libs can be found -----
# ---- choose between lib and lib64 depending on your system ----------
FFTWDIR=$myprefix/fftw-$FFTW_VERSION
FFTWINC=$FFTWDIR/include
FFTWLIB=$FFTWDIR/lib64

# ---- set PFFT install path ---------------------------------------
INSTDIR=$myprefix/pfft-$PFFT_VERSION

# ---- set build directory and name of the log file ------------------
BUILDDIR="/LOCAL/builds/pfft-$PFFT_VERSION"
LOGFILE="$BUILDDIR/build.log"


# ---- set MPI compilers and compiler flags ------------------------
# ---- choose between debugging and optimization flags -------------
COMP="CC=mpicc FC=mpif90 MPICC=mpicc MPIFC=mpif90"
echo "Use PFFT debugging flags? (y/n)"
read answer
if [ ${answer} = "y" ]; then
  INSTDIR="$INSTDIR-dbg"
  BUILDDIR="$BUILDDIR-dbg"
  CFLAGS="-O0 -ggdb -Wall"
  FCFLAGS="-O0 -ggdb -Wall"
  PFFTDBG="--enable-debug"
else
  CFLAGS="-O3 -ffast-math -Wall"
  FCFLAGS="-O3 -Wall"
  PFFTDBG=""
fi

# ---- bash check if directory exists -----------------------------
if [ -d $BUILDDIR ]; then
  echo "Directory $BUILDDIR already exists. Delete it? (y/n)"
  read answer
  if [ ${answer} = "y" ]; then
    rm -rf $BUILDDIR
  else
    echo "Install script aborted."
    exit 1
  fi
fi
mkdir -p $BUILDDIR 

# ---- Normally, we should not have to rerun bootstrap -----------
echo "Rerun bootstrap? (y/n)"
read answer
if [ ${answer} = "y" ]; then
  cd $CONFDIR
  ./bootstrap.sh 2>&1 | tee $LOGFILE
  cd -
fi

# ---- configure, build, install ----------------------------------
cd $BUILDDIR
$CONFDIR/configure --prefix=$INSTDIR --disable-shared \
  $PFFTDBG \
  CPPFLAGS="-I$FFTWINC" LDFLAGS="-L$FFTWLIB" \
  $COMP \
  CFLAGS="$CFLAGS" FCFLAGS="$FCFLAGS" \
  2>&1 | tee -a $LOGFILE
echo "*** PFFT successfully build in"
echo "  $BUILDDIR"

make -j4 check 2>&1 | tee -a $LOGFILE 
make install 2>&1 | tee -a $LOGFILE 

# ---- status output ---------------------------------------------
echo "*** PFFT build in"
echo "  $BUILDDIR"
echo "*** PFFT installed in"
echo "  $INSTDIR"
echo "*** For details, look at the log file"
echo "  $LOGFILE"
