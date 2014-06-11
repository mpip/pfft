#
# Copyright (c) 2008 Jens Keiner
# Copyright (c) 2010,2011 Michael Pippig
# Copyright (c) 2012 Michael Hofmann
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.


AC_DEFUN_ONCE([AX_LIB_PNFFT],[
  _AX_LIB_PNFFT_EXT([],[])
])

AC_DEFUN_ONCE([AX_LIB_PNFFTF],[
  _AX_LIB_PNFFT_EXT([f],[])
])

AC_DEFUN_ONCE([AX_LIB_PNFFTL],[
  _AX_LIB_PNFFT_EXT([l],[])
])


AC_DEFUN_ONCE([AX_LIB_PNFFT_EXT],[
  _AX_LIB_PNFFT_EXT([],[$1])
])

AC_DEFUN_ONCE([AX_LIB_PNFFTF_EXT],[
  _AX_LIB_PNFFT_EXT([f],[$1])
])

AC_DEFUN_ONCE([AX_LIB_PNFFTL_EXT],[
  _AX_LIB_PNFFT_EXT([l],[$1])
])


AC_DEFUN([_AX_LIB_PNFFT_EXT],[

#  AC_REQUIRE([AX_LIB_PFFT]) 

  _AX_LIB_PNFFT_ARGS

  if test "x$2" = "xnobuilt" && test "x${ax_lib_pnfft_with_pnfft}" = "xbuilt" ; then
    ax_lib_pnfft_with_pnfft=
    ax_lib_pnfft_with_pnfft_prefix=
    ax_lib_pnfft_with_pnfft_inc_dir=
    ax_lib_pnfft_with_pnfft_lib_dir=
  fi

  _AX_LIB_PNFFT_CHECK([$1])
])


AC_DEFUN([_AX_LIB_PNFFT_ARGS],[

  AC_ARG_WITH(pnfft,
    [AC_HELP_STRING([--with-pnfft=DIR], [compile with pnfft in DIR])],
    [ax_lib_pnfft_with_pnfft="$withval"], [ax_lib_pnfft_with_pnfft=])

  AC_ARG_WITH(pnfft-prefix,
    [AC_HELP_STRING([--with-pnfft-prefix=PREFIX], [compile with pnfft prefix PREFIX])],
    [ax_lib_pnfft_with_pnfft_prefix="$withval"], [ax_lib_pnfft_with_pnfft_prefix=])

  AC_ARG_WITH(pnfft-includedir,
    [AC_HELP_STRING([--with-pnfft-includedir=DIR], [compile with pnfft include directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_pnfft_with_pnfft_inc_dir="$withval"], [ax_lib_pnfft_with_pnfft_inc_dir=])

  AC_ARG_WITH(pnfft-libdir,
    [AC_HELP_STRING([--with-pnfft-libdir=DIR], [compile with pnfft library directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_pnfft_with_pnfft_lib_dir="$withval"], [ax_lib_pnfft_with_pnfft_lib_dir=])
])


AC_DEFUN([AX_LIB_PNFFT_CHECK],[
  _AX_LIB_PNFFT_CHECK([],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_PNFFTF_CHECK],[
  _AX_LIB_PNFFT_CHECK([f],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_PNFFTL_CHECK],[
  _AX_LIB_PNFFT_CHECK([l],[$1],[$2],[$3],[$4])
])


AC_DEFUN([_AX_LIB_PNFFT_CHECK],[

  ax_type_suffix="$1"
  ax_with_pnfft=m4_ifnblank($2,"$2","$ax_lib_pnfft_with_pnfft")
  ax_with_pnfft_prefix=m4_ifnblank($3,"$3","$ax_lib_pnfft_with_pnfft_prefix")
  ax_with_pnfft_inc_dir=m4_ifnblank($4,"$4","$ax_lib_pnfft_with_pnfft_inc_dir")
  ax_with_pnfft_lib_dir=m4_ifnblank($5,"$5","$ax_lib_pnfft_with_pnfft_lib_dir")

#  echo "with_pnfft: $ax_with_pnfft"
#  echo "with_pnfft_prefix: $ax_with_pnfft_prefix"
#  echo "with_pnfft_lib_dir: $ax_with_pnfft_lib_dir"
#  echo "with_pnfft_inc_dir: $ax_with_pnfft_inc_dir"

  pnfft_CPPFLAGS=
  pnfft_LDFLAGS=
  pnfft_PREFIX=
  pnfft_LIBS=
  pnfft_threads_LIBS=
  pnfft_mpi_LIBS=

  if test "x$ax_with_pnfft" != x ; then
    if test "x$ax_with_pnfft_inc_dir" = x ; then
      ax_with_pnfft_inc_dir="$ax_with_pnfft/include"
    fi
    if test "x$ax_with_pnfft_lib_dir" = x ; then 
      ax_with_pnfft_lib_dir="$ax_with_pnfft/lib"
    fi
  fi

  if test "x$ax_with_pnfft_inc_dir" != x ; then
    for ax_d in $ax_with_pnfft_inc_dir ; do
      if test "x$ax_with_pnfft" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The include directory '$ax_d' does not exist.])])
      fi
      pnfft_CPPFLAGS="$pnfft_CPPFLAGS -I$ax_d"
    done
    pnfft_CPPFLAGS="${pnfft_CPPFLAGS# }"
  fi

  if test "x$ax_with_pnfft_lib_dir" != x ; then 
    for ax_d in $ax_with_pnfft_lib_dir ; do
      if test "x$ax_with_pnfft" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The library directory '$ax_d' does not exist.])])
      fi
      pnfft_LDFLAGS="$pnfft_LDFLAGS -L$ax_d"
    done
    pnfft_LDFLAGS="${pnfft_LDFLAGS# }"
  fi

  AC_LANG_PUSH([C])

  saved_CFLAGS="$CFLAGS"
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS $pfft_CPPFLAGS $pnfft_CPPFLAGS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS $pfft_LDFLAGS $pnfft_LDFLAGS"
  LIBS="$fftw3_mpi_LIBS $fftw3_LIBS $pfft_LIBS $LIBS"

  ax_lib_pnfft=no

  pnfft_PREFIX="$ax_with_pnfft_prefix"

  if test "x$ax_with_pnfft" = xbuilt ; then
#    AC_MSG_NOTICE([pnfft is built!])
    if test -z "$ax_with_pnfft_inc_dir" || test -z "$ax_with_pnfft_lib_dir"; then
      AC_MSG_ERROR([for an in-tree pnfft, set include and library directory with --with-pnfft-include-dir and --with-pnfft-lib-dir])
    fi
    ax_lib_pnfft=yes
    pnfft_LIBS="-l${ax_with_pnfft_prefix}pnfft${ax_type_suffix}"
  else
    # Check if header is present and usable.
    AC_CHECK_HEADER([pnfft.h], [ax_lib_pnfft=yes])

    if test "x$ax_lib_pnfft" = xyes ; then
      saved_LIBS="$LIBS"

      AC_CHECK_LIB([${ax_with_pnfft_prefix}pnfft${ax_type_suffix}], [${ax_with_pnfft_prefix}pnfft${ax_type_suffix}_trafo], [], [ax_lib_pnfft=no])
      test "x${ax_lib_pnfft}" != xno && pnfft_LIBS="-l${ax_with_pnfft_prefix}pnfft${ax_type_suffix}"

      LIBS="$saved_LIBS"
    fi
  fi

  # Restore saved flags.
  CFLAGS="$saved_CFLAGS"
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  AC_LANG_POP([C])

  AC_SUBST(pnfft_CPPFLAGS)
  AC_SUBST(pnfft_LDFLAGS)
  AC_SUBST(pnfft_PREFIX)
  AC_SUBST(pnfft_LIBS)

  if test "x$ax_lib_pnfft" = xyes ; then
    AC_DEFINE_UNQUOTED([PNFFT_PREFIX],[$pnfft_PREFIX],[Define to the prefix of the namespace of the PNFFT library.])
  fi
])
