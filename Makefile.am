## Process this file with automake to create Makefile.in

# Makefile for PFFT.

LIBTOOL_DEPS = @LIBTOOL_DEPS@

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Traverse these subdirectories, the current one first (for config.h).
SUBDIRS = kernel gcell util api .

if ENABLE_TESTS
SUBDIRS += tests
endif

if ENABLE_DOC
SUBDIRS += doc
endif


EXTRA_DIST = bootstrap.sh CONVENTIONS pfft.pc.in

LIBPFFT_OMP_LA =
if ENABLE_OPENMP
LIBPFFT_OMP_LA += lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp.la
endif

# Libraries that are built and installed.
# set library name according to precision
if ENABLE_LIBRARY_INSTALL
lib_LTLIBRARIES = lib@PFFT_PREFIX@pfft@PREC_SUFFIX@.la $(LIBPFFT_OMP_LA)
else
noinst_LTLIBRARIES = lib@PFFT_PREFIX@pfft@PREC_SUFFIX@.la $(LIBPFFT_OMP_LA)
endif


# Libtool convenience libraries.
lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_la_SOURCES =
lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_la_LIBADD = \
	kernel/libkernel.la \
	util/libutil.la \
	api/libapi.la \
	gcell/libgcell.la

lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@


if ENABLE_OPENMP
# Libtool convenience libraries.
lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp_la_SOURCES =
lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp_la_LIBADD = \
	kernel/libkernel.la \
	util/libutil.la \
	api/libapi_omp.la \
	gcell/libgcell.la

lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@
lib@PFFT_PREFIX@pfft@PREC_SUFFIX@_omp_la_CFLAGS = $(OPENMP_CFLAGS)
endif


# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

pfft1@PREC_SUFFIX@.pc: pfft.pc
	cp -f pfft.pc pfft1@PREC_SUFFIX@.pc


#################################################################
# Documentation
#################################################################
if ENABLE_DOC
.PHONY: FORCE doc manual
doc manual: FORCE
	cd doc; $(MAKE) --print-directory $@
FORCE:
endif
