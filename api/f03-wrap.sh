#! /bin/sh

# Script to generate Fortran 2003 wrappers for FFTW's MPI functions.  This
# is necessary because MPI provides no way to deal with C MPI_Comm handles
# from Fortran (where MPI_Comm == integer), but does provide a way to
# deal with Fortran MPI_Comm handles from C (via MPI_Comm_f2c).  So,
# every FFTW function that takes an MPI_Comm argument needs a wrapper
# function that takes a Fortran integer and converts it to MPI_Comm.

# pfft.h depends on fftw3-mpi.h and fftw3.h
# set these paths such that the preprocessor can find the required headers
FFTW_INC=$HOME/local/fftw-3.3.3/include

echo "/* Generated automatically.  DO NOT EDIT! */"
echo

echo "#include \"pfft.h\""
echo "#include \"ipfft.h\""
echo

# Declare prototypes using FFTW_EXTERN, important for Windows DLLs
mpicc -E pfft.h -I${FFTW_INC}  |grep 'pfftl_init' |tr ';' '\n' |grep 'MPI_Comm' |grep -v 'printf' |perl genf03-wrap.pl |grep "MPI_Fint" |sed 's/^/PFFT_EXTERN /;s/$/;/'
mpicc -E pfft.h -I${FFTW_INC}  |grep 'pfftl_init' |tr ';' '\n' |grep 'MPI_Comm' |grep -v 'printf' |perl genf03-wrap.pl

