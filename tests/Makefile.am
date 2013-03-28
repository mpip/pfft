# Subdirectories
SUBDIRS = 

# Enable Fortran subdirectory only on request.
if ENABLE_FORTRAN
SUBDIRS += fortran f03
endif

# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@.la $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = \
	minimal_test \
	minimal_check_c2c minimal_check_c2c_transposed

check_PROGRAMS += \
	simple_check_c2c simple_check_c2c_transposed \
	simple_check_r2c simple_check_r2c_transposed \
	simple_check_r2r simple_check_r2r_transposed



#if ENABLE_COMMON_PFFT
# do not compile checks
#else
check_PROGRAMS += \
	simple_check_c2c_4d simple_check_c2c_4d_transposed \
	simple_check_r2c_4d simple_check_r2c_4d_transposed \
	simple_check_r2r_4d simple_check_r2r_4d_transposed

check_PROGRAMS += \
       	simple_check_c2c_4d_on_3d simple_check_c2c_4d_on_3d_transposed \
       	simple_check_r2c_4d_on_3d simple_check_r2c_4d_on_3d_transposed \
       	simple_check_r2r_4d_on_3d simple_check_r2r_4d_on_3d_transposed

check_PROGRAMS += \
       	simple_check_c2c_3d_on_3d simple_check_c2c_3d_on_3d_transposed \
       	simple_check_r2c_3d_on_3d simple_check_r2c_3d_on_3d_transposed \
       	simple_check_r2r_3d_on_3d simple_check_r2r_3d_on_3d_transposed

check_PROGRAMS += \
	simple_check_ousam_c2c simple_check_ousam_c2c_transposed \
	simple_check_ousam_r2c simple_check_ousam_r2c_transposed \
	simple_check_ousam_c2c_4d simple_check_ousam_c2c_4d_transposed

check_PROGRAMS += \
	simple_check_ousam_r2c_4d simple_check_ousam_r2c_4d_transposed \
	simple_check_ousam_c2c_4d_on_3d simple_check_ousam_c2c_4d_on_3d_transposed \
	simple_check_ousam_r2c_4d_on_3d simple_check_ousam_r2c_4d_on_3d_transposed

check_PROGRAMS += \
	simple_check_ghost_c2c \
	simple_check_ghost_c2c_3d_on_3d 

check_PROGRAMS += \
	time_c2c_3d_on_3d_transposed \
	time_c2c_transposed
#endif

