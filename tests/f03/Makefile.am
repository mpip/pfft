# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Directory of pnfft.f (which is build first)
AM_FCCPPFLAGS = -I$(top_builddir)/api

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@.la $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# noinst_LTLIBRARIES = libtests.la

# search for Fortran source files to build objects
AM_DEFAULT_SOURCE_EXT = .f90

check_PROGRAMS = 

# These programs are built by 'make check' and may be tested afterwards.
#if ENABLE_COMMON_PNFFT
# do not compile checks
#else
check_PROGRAMS += \
	simple_check_c2c \
	simple_check_ousam_c2c_4d_on_3d
#endif

