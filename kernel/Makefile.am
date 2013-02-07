# Directory of ipfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel

# Directory of pfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# Directory of util.h
AM_CPPFLAGS += -I$(top_srcdir)/util

noinst_LTLIBRARIES = libkernel.la

# Group local sources into convenience lib.
libkernel_la_SOURCES = \
	block.c \
	check.c \
	ipfft.h \
	malloc.c \
	ousample.c \
	outrafo.c \
	partrafo-transposed.c \
	partrafo.c \
	sertrafo.c \
	procmesh.c \
	remap_3dto2d.c \
	timer.c \
	transpose.c
