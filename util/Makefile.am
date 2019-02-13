# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api 
#
# Directory of fftw3
AM_CPPFLAGS += -I$(top_srcdir)/fftw3/api
AM_CPPFLAGS += -I$(top_srcdir)/fftw3/mpi

noinst_LTLIBRARIES = libutil.la

# Group local sources into convenience lib.
libutil_la_SOURCES = \
	util.c \
	util.h \
	getargs.c

