# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Directory of pnfft.f (which is build first)
AM_FCCPPFLAGS = -I$(top_builddir)/api
#
# Directory of fftw3
AM_FCCPPFLAGS += -I$(top_srcdir)/fftw3/api
AM_FCCPPFLAGS += -I$(top_srcdir)/fftw3/mpi

# OpenMP support
AM_FCFLAGS = $(OPENMP_FCFLAGS)

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@@OPENMP_SUFFIX@.la $(fftw3_openmp_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)
LDADD += $(top_builddir)/fftw3/libfftw3@PREC_SUFFIX@@OPENMP_SUFFIX@.la
LDADD += $(top_builddir)/fftw3/mpi/libfftw3@PREC_SUFFIX@_mpi.la

# noinst_LTLIBRARIES = libtests.la

# search for Fortran source files to build objects
AM_DEFAULT_SOURCE_EXT = .F03

check_PROGRAMS = 

# These programs are built by 'make check' and may be tested afterwards.
#if ENABLE_COMMON_PNFFT
# do not compile checks
#else
check_PROGRAMS += \
	simple_check_c2c \
	simple_check_c2c_transposed \
	simple_check_r2c_4d_on_3d \
	simple_check_r2c_4d_on_3d_transposed \
	simple_check_ousam_c2c_4d_on_3d \
	simple_check_ousam_c2c_4d_on_3d_transposed \
	simple_check_ousam_r2c_4d_on_3d \
	simple_check_ousam_r2c_4d_on_3d_transposed \
	time_c2c \
	time_c2c_transposed
#endif

