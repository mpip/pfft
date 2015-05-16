# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# OpenMP support
AM_CFLAGS = $(OPENMP_CFLAGS)

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp.la $(fftw3_openmp_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = \
	simple_check_c2c simple_check_c2c_transposed\
	simple_check_c2c_omp omp_bench_c2c\
	transpose_check_fftw
