# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Directory of pfft.f (which is build first)
AM_FCCPPFLAGS = -I$(top_builddir)/api

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@.la $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# noinst_LTLIBRARIES = libtests.la

# search for Fortran source files to build objects
AM_DEFAULT_SOURCE_EXT = .F90

# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = \
	simple_test

#if ENABLE_COMMON_PFFT
# do not compile checks
#else
check_PROGRAMS += \
	simple_check_c2c simple_check_c2c_transposed \
	simple_check_r2c simple_check_r2c_transposed \
	simple_check_r2r simple_check_r2r_transposed

check_PROGRAMS += \
	simple_check_r2c_3d_on_3d simple_check_r2c_3d_on_3d_transposed

check_PROGRAMS += \
	simple_check_ghost_c2c simple_check_ghost_c2c_3d_on_3d
#endif

