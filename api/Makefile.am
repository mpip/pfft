# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of util.h
AM_CPPFLAGS += -I$(top_srcdir)/util

# Headers that are installed.
if ENABLE_HEADER_INSTALL
include_HEADERS = pfft.h
nodist_include_HEADERS = pfft.f
else
noinst_HEADERS = pfft.h
endif
	
noinst_LTLIBRARIES = libapi.la

# Group local sources into convenience lib.
libapi_la_SOURCES = \
        api-basic.c \
        api-adv.c \
	pfft.f.in \
	fortran-api.c \
	fortran-mangling.h \
	fortran-prototypes.h \
	fortran-wrappers.h

