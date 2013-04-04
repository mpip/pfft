#
# Copyright (c) 2011-2013 Michael Pippig
# Copyright (c) 2005, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
#
# This file is part of PFFT.
#
# PFFT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PFFT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PFFT.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# Parts of this file are based on the configure.ac script from NFFT 3.2.0-alpha
# and the configure.ac script from FFTW 3.3.1-alpha 
################################################################################
# Process this file with autoconf to produce a configure script.
################################################################################
# NOTE: If you just want to build PFFT, do not use this file. Just follow the
# installation instructions as described in the user manual
################################################################################

# Initialize Autoconf. We require version 2.59 or newer.
AC_PREREQ(2.64)

# autoconf initialization
AC_INIT([PFFT],[1.0.7-alpha],[michael.pippig@mathematik.tu-chemnitz.de],[pfft],[http://www.tu-chemnitz.de/~mpip/software/])
m4_ifndef([AC_PACKAGE_URL],
	  [AC_SUBST([PACKAGE_URL], [http://www.tu-chemnitz.de/~mpip/software/])])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring in common/pfft                              *])
AC_MSG_NOTICE([****************************************************************])

# copyright notice
AC_COPYRIGHT([Copyright (C) 2011-2013 Michael Pippig])

# m4 macros go here
AC_CONFIG_MACRO_DIR([m4])

# where to put auxilliary files
AC_CONFIG_AUX_DIR([build-aux])

# canonical host system type string
AC_CANONICAL_HOST

# automake initialization; require 1.9.6 or newer, do not enforce GNU rules.
AM_INIT_AUTOMAKE([1.9.6 foreign])

################################################################################
# options for customizing the build process
################################################################################

# multithreaded code
AC_ARG_ENABLE(threads, [AS_HELP_STRING([--enable-threads],
  [enable multithreaded code])], enable_threads=$enableval, enable_threads=no)
AM_CONDITIONAL(HAVE_THREADS, test "x$enable_threads" = "xyes" )

# debug mode
AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],
  [compile with extra runtime checks for debugging])], enable_debug=$enableval,
  enable_debug=no)
if test "x$enable_debug" = "xyes"; then
  AC_DEFINE(PFFT_DEBUG,1,[Define to enable extra debugging code.])
fi

# disable fortran wrappers
AC_ARG_ENABLE(fortran, [AS_HELP_STRING([--disable-fortran],
  [do not include Fortran-callable wrappers])], enable_fortran=$enableval,
  enable_fortran=yes)

# include gfortran wrappers (although we use an other compiler)
AC_ARG_WITH(gfortran-wrappers, [AS_HELP_STRING([--with-gfortran-wrappers],
  [force inclusion of gfortran-compatible wrappers in addition to any other Fortran compiler that is detected])],
  with_gfortran_wrappers=$withval)

# set precision
AC_ARG_ENABLE(single, [AS_HELP_STRING([--enable-single],[compile pfft in single precision])], ok=$enableval, ok=no)
AC_ARG_ENABLE(float,  [AS_HELP_STRING([--enable-float], [synonym for --enable-single])], ok=$enableval)
if test "$ok" = "yes"; then
  AC_DEFINE(PFFT_SINGLE,1,[Define to compile in single precision.])
  PRECISION=s
else
  PRECISION=d
fi
AM_CONDITIONAL(SINGLE, test "$ok" = "yes")

AC_ARG_ENABLE(long-double, [AS_HELP_STRING([--enable-long-double],[compile pfft in long-double precision])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  if test "$PRECISION" = "s"; then
    AC_MSG_ERROR([--enable-single/--enable-long-double conflict])
  fi
  AC_DEFINE(PFFT_LDOUBLE,1,[Define to compile in long-double precision.])
  PRECISION=l
fi
AM_CONDITIONAL(LDOUBLE, test "$ok" = "yes")
AC_SUBST(PRECISION)

# compute library suffix
case "$PRECISION" in
     s) PREC_SUFFIX=f;;
     d) PREC_SUFFIX=;;
     l) PREC_SUFFIX=l;;
esac
AC_SUBST(PREC_SUFFIX)

# Add options (--disable-header-install, --disable-library-install, and --disable-tests) and Automake conditionals to disable install of headers and libraries and build of tests on request.
# If the AX_DISABLE macro is missing, then headers, libraries, and tests will be enabled on default.
m4_ifdef([AX_DISABLE],[AX_DISABLE([header-install,library-install,tests])])
AM_CONDITIONAL([ENABLE_HEADER_INSTALL], [test "x$ax_disable_header_install" != xyes])
AM_CONDITIONAL([ENABLE_LIBRARY_INSTALL], [test "x$ax_disable_library_install" != xyes])
AM_CONDITIONAL([ENABLE_TESTS], [test "x$ax_disable_tests" != xyes])

# library prefix
AC_ARG_ENABLE(lib-prefix, [AS_HELP_STRING([--enable-lib-prefix=PREFIX],[use PREFIX for naming the pfft library])], PFFT_PREFIX=$enableval, PFFT_PREFIX=)
AC_SUBST(PFFT_PREFIX)


################################################################################
# C compiler characteristis
################################################################################

# select programming language
AC_LANG(C)

AX_PROG_CC_MPI(,,AC_MSG_FAILURE([PFFT requires an MPI C compiler.]))

# The AC_PROG_CC_C99 macro is relatively new in Autoconf, so
# if the developer only used an older version, the user will
# have to specify the needed flags (e.g., CC="gcc -std=gnu99").
m4_ifdef([AC_PROG_CC_C99], [AC_PROG_CC_C99])

# enable "inline" keyword
AC_C_INLINE


################################################################################
# Fortran compiler characteristics
################################################################################

# set temporary programming language to Fortran
AC_LANG_PUSH([Fortran])

if test "$enable_fortran" = "yes"; then
  # Get the Fortran MPI compiler.
  AX_PROG_FC_MPI(,,[
		AC_MSG_WARN([*** Could not find MPI Fortran compiler; using default Fortran wrappers.])
		enable_fortran=no
  ])
else
  AC_DEFINE([DISABLE_FORTRAN], 1, [Define to disable Fortran wrappers.])
fi

AM_CONDITIONAL(ENABLE_FORTRAN,[test "x$enable_fortran" = "xyes"])

if test "x$enable_fortran" = xyes; then
   # Find out how to name a MAIN Fortran entry point in C.
   AC_FC_MAIN([], [enable_fortran=no
        AC_MSG_WARN([*** Could not figure out how to link C and Fortran; using default Fortran wrappers.])])

  # Get it to accept long source file lines.
  AC_FC_LINE_LENGTH
  # Find out how to link against Fortran libraries.
  AC_FC_LIBRARY_LDFLAGS
  # Bugfix for BlueGene xlf2003: remove spurious -link entry.
  ACX_FC_LIBRARY_LDFLAGS_FIX
  # Get the Fortran compiler to accept free-form source.
  AC_FC_FREEFORM
  # Accept preprocessed source files names ending in ".F90".
  AC_FC_PP_SRCEXT([F90])
  # Get the flag to use for setting Fortran defines on the command line.
  AC_FC_PP_DEFINE
  # Find out how to add to the Fortran include path.
  AC_FC_MODULE_FLAG
  # Find out the Fortran module file extension.
  AC_FC_MODULE_EXTENSION

  # Set FCDEFS and adjust -D flag if needed.
  FCDEFS="$DEFS -DHAVE_FCONFIG_H"
  if test "x$FC_DEFINE" != x-D; then
    FCDEFS=`echo " $FCDEFS" | sed "s@ -D@ $FC_DEFINE@g"`
  fi
  AC_SUBST([FCDEFS])

  # Find out how to call Fortran functions from C.
  AC_FC_WRAPPERS
  
  AC_FC_FUNC(fcfoo)
  AC_FC_FUNC(fc_foo)
  fc_foo2=`echo $fcfoo | sed 's/c/c_/'`
  if test "$fc_foo" = "$fc_foo2"; then
    AC_DEFINE(FC_FUNC_EQUIV, 1, [Define if FC_FUNC and FC_FUNC_ are equivalent.])

    # Include g77 wrappers by default for GNU systems or gfortran
    with_gfortran_wrappers=$ac_cv_fc_compiler_gnu
    case $host_os in *gnu*) with_gfortran_wrappers=yes ;; esac
  fi
else
  with_gfortran_wrappers=no
fi

if test "x$with_gfortran_wrappers" = "xyes"; then
  AC_DEFINE(WITH_GFORTRAN_WRAPPERS,1,[Include gfortran-compatible wrappers in addition to any other Fortran wrappers.])
fi

# return to former programming language
AC_LANG_POP([Fortran])


################################################################################
# Libtool initialization
################################################################################

# minimum required libtool version
LT_PREREQ([2.2.0])

# libtool initialization
LT_INIT([win32-dll])

# substitute LIBTOOL_DEPS variable in generated Makefiles
AC_SUBST([LIBTOOL_DEPS])

# version information for shared library
SHARED_VERSION_INFO="0:0:0"

# substitute SHARED_VERSION_INFO in generated Makefiles
AC_SUBST(SHARED_VERSION_INFO)


################################################################################
# program checks
################################################################################

# check to build a library
# Replaced by LT_INIT
# AC_PROG_RANLIB


################################################################################
# 3rd party libraries
################################################################################

# May need sincos from libm.
AC_CHECK_LIB([m], [sincos])

# Check for OpenMP.
AX_OPENMP

# Check for fftw3, MPI FFTW and threaded FFTW.
if test "x$PRECISION" = "xs" ; then
  AX_LIB_FFTW3F
elif test "x$PRECISION" = "xl" ; then
  AX_LIB_FFTW3L
else
  AX_LIB_FFTW3
fi

if test "x$ax_lib_fftw3" = "xno"; then
  AC_MSG_ERROR([You do not seem to have the FFTW-3.3 library installed. You can ]
      [download it from http://www.fftw.org. If you have installed FFTW-3.3, ]
      [make sure that this configure script can find it. See ./configure ]
      [--help for more information.])
fi

if test "x$ax_lib_fftw3_mpi" = "xno"; then
  AC_MSG_ERROR([You do not seem to have the MPI part of the  FFTW-3.3 library installed.]
      [You can download it from http://www.fftw.org. If you have installed FFTW-3.3, ]
      [make sure that this configure script can find it. See ./configure --help]
      [for more information.])
fi

if  test "x$enable_threads" = "xyes" -a "x$ax_lib_fftw3_threads" = "xno"; then
  AC_MSG_ERROR([You do not seem to have the threaded FFTW-3.3 library installed.])
fi


################################################################################
# compiler options
################################################################################

# Try to choose good compiler options.
if test "x$ac_test_CFLAGS" != "xset"; then
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"
  LIBS="$LIBS $fftw3_LIBS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS"
  AX_CC_MAXOPT
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"
fi

CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"
LDFLAGS="$LDFLAGS $fftw3_LDFLAGS"

# It is probably confusing for the user to have both CPPFLAGS and FCCPPFLAGS
# (but it may be necessary to set them differently in some situations).  Since
# we don't use any preprocessor tests below, default the latter to the former
# at 'make' run time if they are equal, so setting the former affects both.
if test "x$CPPFLAGS" = "x$FCCPPFLAGS" || test -z "$FCCPPFLAGS"; then
  FCCPPFLAGS=`echo " $CPPFLAGS" | sed "s@ -D@ $FC_DEFINE@g"`
fi
AC_SUBST([FCCPPFLAGS])

# add debug flags to CFLAGS and FCFLAGS, in debug mode only
if test "x${enable_debug}" = "xyes";
then
  # add debug symbols for all compilers
  CFLAGS="$CFLAGS -g"
  FCFLAGS="$FCFLAGS -g"
  
  # add gcc warnings
  if test "x$ac_cv_c_compiler_gnu" = "xyes"; then
    CFLAGS="$CFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
    CFLAGS="$CFLAGS -Wno-long-long -Wshadow -Wbad-function-cast -Wwrite-strings"
    CFLAGS="$CFLAGS -Wstrict-prototypes -Wredundant-decls -Wnested-externs"
    CFLAGS="$CFLAGS -Wundef -Wconversion -Wmissing-prototypes "
    CFLAGS="$CFLAGS -Wmissing-declarations"
    FCFLAGS="$FCFLAGS -Wall -W -Wcast-align -pedantic"
    FCFLAGS="$FCFLAGS -Wshadow -Wconversion"
    # MPI defines a lot of parameters, suppress the warnings
    FCFLAGS="$FCFLAGS -Wno-unused-parameter"
    # PFFT uses large integer for 64-bit support, 
    # turn off annoying conversion warnings
    CFLAGS="$CFLAGS -Wno-conversion"
    FCFLAGS="$FCFLAGS -Wno-conversion"
  fi
  if test "x$ac_cv_c_compiler_ibm" = "xyes"; then
    CFLAGS="$CFLAGS -qfullpath -qcheck=all"
    FCFLAGS="$FCFLAGS -qfullpath"
  fi
fi

# option to accept C99
CFLAGS="$CFLAGS $ac_cv_prog_cc_c99" 


################################################################################
# header files/data types/compiler characteristics
################################################################################

# AC_CHECK_HEADERS([mpi.h stdarg.h])

# AC_CHECK_HEADERS([math.h stdio.h stdlib.h time.h  sys/time.h \
#   complex.h string.h float.h limits.h stdarg.h stddef.h sys/types.h stdint.h \
#   inttypes.h stdbool.h malloc.h c_asm.h intrinsics.h mach/mach_time.h])

# ToDo: pfft depends on ptrdiff_t only in the header, define INT macro to internal big integer type
# AC_CHECK_TYPES([size_t], [], [AC_MSG_ERROR([size_t is not a supported type with your compiler.])])
AC_CHECK_TYPES([ptrdiff_t], [], [AC_MSG_ERROR([ptrdiff_t is not a supported type with your compiler.])])
AC_CHECK_TYPE([long double],
              [AC_DEFINE(HAVE_LONG_DOUBLE, 1, [Define to 1 if the compiler supports 'long double'])],
[
if test $PRECISION = l; then
  AC_MSG_ERROR([long double is not a supported type with your compiler.])
fi
])

AC_CHECK_SIZEOF([ptrdiff_t])
AC_SUBST([PTRDIFF_T_KIND], ["$ac_cv_sizeof_ptrdiff_t"])

AC_CHECK_SIZEOF(MPI_Fint, [], [#include <mpi.h>])
if test 0 = $ac_cv_sizeof_MPI_Fint; then
   AC_MSG_WARN([sizeof(MPI_Fint) test failed]);
   dnl As a backup, assume Fortran integer == C int
   AC_CHECK_SIZEOF(int) 
   if test 0 = $ac_cv_sizeof_int; then AC_MSG_ERROR([sizeof(int) test failed]); fi
   ac_cv_sizeof_MPI_Fint=$ac_cv_sizeof_int
fi
C_MPI_FINT=C_INT`expr $ac_cv_sizeof_MPI_Fint \* 8`_T
AC_SUBST(C_MPI_FINT)

# We need the math library.
AC_SEARCH_LIBS([log], [m])

# pass location of FFTW to make distcheck
# TODO: this is not a good way to pass the FFTW locations, since it sets
# CFLAGS and LDFLAGS. Try to use --with-fftw3-includedir=... and
# --with-fftw3-libdir=... instead!
DISTCHECK_CONFIGURE_FLAGS="CPPFLAGS=\"$CPPFLAGS\" LDFLAGS=\"$LDFLAGS\""

AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

################################################################################
# create output files
################################################################################

# Create config headers: config.h comes from autoheader-generated config.h.in,
# fconfig.h.in is written manually and only contains stuff needed for Fortran.
AC_CONFIG_HEADERS([config.h])

# Convert config.h into an fconfig.h suitable for inclusion into Fortran.
# The sed script transforms C style /* comments */ into Fortran style.
AC_CONFIG_HEADERS([fconfig.h],
 [[sed '
    /\/\*\(.*\)\*\//{
      s,,!\1,p
      d
    }
    /\/\*/,/\*\//{
      s,/\*,,
      s,\*/,,
      s,^,!,p
      d
    }' fconfig.h > fconfig.ht
   mv -f fconfig.ht fconfig.h
]])

# The output files to be generated
AC_CONFIG_FILES(Makefile \
                pfft.pc \
                api/pfft.f \
                api/pfft.f03 \
                api/pfftl.f03 \
		kernel/Makefile \
		util/Makefile \
		api/Makefile \
		gcell/Makefile \
		tests/Makefile \
		tests/fortran/Makefile \
		tests/f03/Makefile)


AC_OUTPUT
