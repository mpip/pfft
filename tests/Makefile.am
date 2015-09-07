# Subdirectories
SUBDIRS = 

# Enable Fortran subdirectory only on request and only for double precision.
if ENABLE_FORTRAN
if DOUBLE
SUBDIRS += fortran f03
endif
endif

if ENABLE_OPENMP
if DOUBLE
SUBDIRS += openmp
endif
endif

# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# OpenMP support
AM_CFLAGS = $(OPENMP_CFLAGS)

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PFFT_PREFIX@pfft@PREC_SUFFIX@@OPENMP_SUFFIX@.la $(fftw3_openmp_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)


EXTRA_DIST = \
	build_checks.sh \
	run_checks.sh

%.c: %.ctex
	cat $< \
	    | sed 's/\\red§\([^§]*\)§/\1/g; s/\\del§\([^§]*\)§/\1/g # erase single line \red and \del commands' \
	    | sed '/^\\red§$$/d; /^\\del§$$/d; /^§$$/d              # erase multi-line \red and \del commands' \
	    > $@

manual_FILES = man_fftw_3don1d.ctex

EXTRA_DIST += $(manual_FILES)



# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = 

if SINGLE
check_PROGRAMS += \
	simple_check_c2c_float
endif

if LDOUBLE
check_PROGRAMS += \
	simple_check_c2c_ldouble
endif

if DOUBLE
check_PROGRAMS += \
	man_fftw_3don1d

# nodist_man_fftw_3don1d_SOURCES = man_fftw_3don1d.c

check_PROGRAMS += \
	manual_c2c_3d \
	minimal_check_c2c minimal_check_c2c_transposed

check_PROGRAMS += \
	simple_check_c2c_inplace simple_check_c2c_transposed_inplace

check_PROGRAMS += \
	simple_check_c2c simple_check_c2c_transposed \
	simple_check_c2r_c2c simple_check_c2r_c2c_shifted simple_check_c2r_c2c_ousam_shifted \
	simple_check_r2c simple_check_r2c_transposed \
	simple_check_r2r simple_check_r2r_transposed

check_PROGRAMS += \
	simple_check_c2c_newarray simple_check_c2c_transposed_newarray \
	simple_check_r2c_newarray simple_check_r2c_transposed_newarray \
	simple_check_r2c_padded_newarray simple_check_r2c_padded_transposed_newarray \
	simple_check_ousam_c2c_4d_newarray simple_check_ousam_c2c_4d_transposed_newarray \
	simple_check_ousam_r2c_4d_newarray simple_check_ousam_r2c_4d_transposed_newarray

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
	simple_check_r2c_3d_on_3d_transposed_many

check_PROGRAMS += \
	simple_check_ousam_c2c simple_check_ousam_c2c_transposed \
	simple_check_ousam_c2r \
	simple_check_ousam_r2c simple_check_ousam_r2c_transposed \
	simple_check_ousam_r2r simple_check_ousam_r2r_transposed \
	simple_check_ousam_r2c_padded \
	simple_check_ousam_c2c_4d simple_check_ousam_c2c_4d_transposed

check_PROGRAMS += \
	simple_check_ousam_r2c_4d simple_check_ousam_r2c_4d_transposed \
	simple_check_ousam_c2c_4d_on_3d simple_check_ousam_c2c_4d_on_3d_transposed \
	simple_check_ousam_r2c_4d_on_3d simple_check_ousam_r2c_4d_on_3d_transposed

check_PROGRAMS += \
	simple_check_ghost_c2c \
	simple_check_ghost_c2c_3d_on_3d \
	simple_check_ghost_r2c_input \
	simple_check_ghost_r2c_input_padded \
	simple_check_ghost_r2c_output

check_PROGRAMS += \
	bench_c2c \
	time_c2c \
	time_c2c_transposed
endif

