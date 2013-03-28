# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of util.h
AM_CPPFLAGS += -I$(top_srcdir)/util

# Headers that are installed.
if ENABLE_HEADER_INSTALL
include_HEADERS = pfft.h
# provide separate header for long precision, since some compilers do not support it
nodist_include_HEADERS = pfft.f pfft.f03 pfftl.f03
else
noinst_HEADERS = pfft.h
endif
	
noinst_LTLIBRARIES = libapi.la

EXTRA_DIST = f03-api.sh f03-wrap.sh genf03-api.pl genf03-wrap.pl build-f03-api.sh

# Group local sources into convenience lib.
libapi_la_SOURCES = \
        api-basic.c \
        api-adv.c \
	pfft.f.in \
	fortran-api.c \
	fortran-mangling.h \
	fortran-prototypes.h \
	fortran-wrappers.h \
	pfft.f03.in \
	pfftl.f03.in \
	f03-wrap.c

