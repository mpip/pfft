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


AC_DEFUN_ONCE([AX_LIB_PFFT],[
  _AX_LIB_PFFT_EXT([],[])
])

AC_DEFUN_ONCE([AX_LIB_PFFTF],[
  _AX_LIB_PFFT_EXT([f],[])
])

AC_DEFUN_ONCE([AX_LIB_PFFTL],[
  _AX_LIB_PFFT_EXT([l],[])
])


AC_DEFUN_ONCE([AX_LIB_PFFT_EXT],[
  _AX_LIB_PFFT_EXT([],[$1])
])

AC_DEFUN_ONCE([AX_LIB_PFFTF_EXT],[
  _AX_LIB_PFFT_EXT([f],[$1])
])

AC_DEFUN_ONCE([AX_LIB_PFFTL_EXT],[
  _AX_LIB_PFFT_EXT([l],[$1])
])


AC_DEFUN([_AX_LIB_PFFT_EXT],[

#  AC_REQUIRE([AX_LIB_FFTW3]) 

  _AX_LIB_PFFT_ARGS

  if test "x$2" = "xnobuilt" && test "x${ax_lib_pfft_with_pfft}" = "xbuilt" ; then
    ax_lib_pfft_with_pfft=
    ax_lib_pfft_with_pfft_prefix=
    ax_lib_pfft_with_pfft_inc_dir=
    ax_lib_pfft_with_pfft_lib_dir=
  fi

  _AX_LIB_PFFT_CHECK([$1])
])


AC_DEFUN([_AX_LIB_PFFT_ARGS],[

  AC_ARG_WITH(pfft,
    [AC_HELP_STRING([--with-pfft=DIR], [compile with pfft in DIR])],
    [ax_lib_pfft_with_pfft="$withval"], [ax_lib_pfft_with_pfft=])

  AC_ARG_WITH(pfft-prefix,
    [AC_HELP_STRING([--with-pfft-prefix=PREFIX], [compile with pfft prefix PREFIX])],
    [ax_lib_pfft_with_pfft_prefix="$withval"], [ax_lib_pfft_with_pfft_prefix=])

  AC_ARG_WITH(pfft-includedir,
    [AC_HELP_STRING([--with-pfft-includedir=DIR], [compile with pfft include directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_pfft_with_pfft_inc_dir="$withval"], [ax_lib_pfft_with_pfft_inc_dir=])

  AC_ARG_WITH(pfft-libdir,
    [AC_HELP_STRING([--with-pfft-libdir=DIR], [compile with pfft library directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_pfft_with_pfft_lib_dir="$withval"], [ax_lib_pfft_with_pfft_lib_dir=])
])


AC_DEFUN([AX_LIB_PFFT_CHECK],[
  _AX_LIB_PFFT_CHECK([],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_PFFTF_CHECK],[
  _AX_LIB_PFFT_CHECK([f],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_PFFTL_CHECK],[
  _AX_LIB_PFFT_CHECK([l],[$1],[$2],[$3],[$4])
])


AC_DEFUN([_AX_LIB_PFFT_CHECK],[

  ax_type_suffix="$1"
  ax_with_pfft=m4_ifnblank($2,"$2","$ax_lib_pfft_with_pfft")
  ax_with_pfft_prefix=m4_ifnblank($3,"$3","$ax_lib_pfft_with_pfft_prefix")
  ax_with_pfft_inc_dir=m4_ifnblank($4,"$4","$ax_lib_pfft_with_pfft_inc_dir")
  ax_with_pfft_lib_dir=m4_ifnblank($5,"$5","$ax_lib_pfft_with_pfft_lib_dir")

#  echo "with_pfft: $ax_with_pfft"
#  echo "with_pfft_prefix: $ax_with_pfft_prefix"
#  echo "with_pfft_lib_dir: $ax_with_pfft_lib_dir"
#  echo "with_pfft_inc_dir: $ax_with_pfft_inc_dir"

  pfft_CPPFLAGS=
  pfft_LDFLAGS=
  pfft_PREFIX=
  pfft_LIBS=
  pfft_threads_LIBS=
  pfft_mpi_LIBS=

  if test "x$ax_with_pfft" != x ; then
    if test "x$ax_with_pfft_inc_dir" = x ; then
      ax_with_pfft_inc_dir="$ax_with_pfft/include"
    fi
    if test "x$ax_with_pfft_lib_dir" = x ; then 
      ax_with_pfft_lib_dir="$ax_with_pfft/lib"
    fi
  fi

  if test "x$ax_with_pfft_inc_dir" != x ; then
    for ax_d in $ax_with_pfft_inc_dir ; do
      if test "x$ax_with_pfft" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The include directory '$ax_d' does not exist.])])
      fi
      pfft_CPPFLAGS="$pfft_CPPFLAGS -I$ax_d"
    done
    pfft_CPPFLAGS="${pfft_CPPFLAGS# }"
  fi

  if test "x$ax_with_pfft_lib_dir" != x ; then 
    for ax_d in $ax_with_pfft_lib_dir ; do
      if test "x$ax_with_pfft" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The library directory '$ax_d' does not exist.])])
      fi
      pfft_LDFLAGS="$pfft_LDFLAGS -L$ax_d"
    done
    pfft_LDFLAGS="${pfft_LDFLAGS# }"
  fi

  AC_LANG_PUSH([C])

  saved_CFLAGS="$CFLAGS"
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"
  CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS $pfft_CPPFLAGS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS $pfft_LDFLAGS"
  LIBS="$fftw3_mpi_LIBS $fftw3_LIBS $LIBS"

  ax_lib_pfft=no

  pfft_PREFIX="$ax_with_pfft_prefix"

  if test "x$ax_with_pfft" = xbuilt ; then
#    AC_MSG_NOTICE([pfft is built!])
    if test -z "$ax_with_pfft_inc_dir" || test -z "$ax_with_pfft_lib_dir"; then
      AC_MSG_ERROR([for an in-tree pfft, set include and library directory with --with-pfft-include-dir and --with-pfft-lib-dir])
    fi
    ax_lib_pfft=yes
    pfft_LIBS="-l${ax_with_pfft_prefix}pfft${ax_type_suffix}"
  else
    # Check if header is present and usable.
    AC_CHECK_HEADER([pfft.h], [ax_lib_pfft=yes])

    if test "x$ax_lib_pfft" = xyes ; then
      saved_LIBS="$LIBS"

      AC_CHECK_LIB([${ax_with_pfft_prefix}pfft${ax_type_suffix}], [${ax_with_pfft_prefix}pfft${ax_type_suffix}_execute], [], [ax_lib_pfft=no])
      test "x${ax_lib_pfft}" != xno && pfft_LIBS="-l${ax_with_pfft_prefix}pfft${ax_type_suffix}"

      LIBS="$saved_LIBS"
    fi
  fi

  # Restore saved flags.
  CFLAGS="$saved_CFLAGS"
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  AC_LANG_POP([C])

  AC_SUBST(pfft_CPPFLAGS)
  AC_SUBST(pfft_LDFLAGS)
  AC_SUBST(pfft_PREFIX)
  AC_SUBST(pfft_LIBS)

  if test "x$ax_lib_pfft" = xyes ; then
    AC_DEFINE_UNQUOTED([PFFT_PREFIX],[$pfft_PREFIX],[Define to the prefix of the namespace of the PFFT library.])
  fi
])
