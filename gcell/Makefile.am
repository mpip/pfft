# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api 

# Directory of util.h
AM_CPPFLAGS += -I$(top_srcdir)/util

noinst_LTLIBRARIES = libgcell.la

# Group local sources into convenience lib.
libgcell_la_SOURCES = \
	gcells_plan.c \
	gcells_RMA.c \
	gcells_sendrecv.c \
	gctimer.c
