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


AC_DEFUN_ONCE([AX_LIB_FFTW3],[
  _AX_LIB_FFTW3_EXT([],[])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3F],[
  _AX_LIB_FFTW3_EXT([f],[])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3L],[
  _AX_LIB_FFTW3_EXT([l],[])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3Q],[
  _AX_LIB_FFTW3_EXT([q],[])
])


AC_DEFUN_ONCE([AX_LIB_FFTW3_EXT],[
  _AX_LIB_FFTW3_EXT([],[$1])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3F_EXT],[
  _AX_LIB_FFTW3_EXT([f],[$1])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3L_EXT],[
  _AX_LIB_FFTW3_EXT([l],[$1])
])

AC_DEFUN_ONCE([AX_LIB_FFTW3Q_EXT],[
  _AX_LIB_FFTW3_EXT([q],[$1])
])


AC_DEFUN([_AX_LIB_FFTW3_EXT],[

#  AC_REQUIRE([AX_OPENMP]) 

  # Check parameters.
  _AX_LIB_FFTW3_ARGS

  if test "x$2" = "xnobuilt" && test "x${ax_lib_fftw3_with_fftw3}" = "xbuilt" ; then
    ax_lib_fftw3_with_fftw3=
    ax_lib_fftw3_with_fftw3_prefix=
    ax_lib_fftw3_with_fftw3_inc_dir=
    ax_lib_fftw3_with_fftw3_lib_dir=
  fi

  _AX_LIB_FFTW3_CHECK([$1])
])


AC_DEFUN([_AX_LIB_FFTW3_ARGS],[

  AC_ARG_WITH(fftw3,
    [AC_HELP_STRING([--with-fftw3=DIR], [compile with fftw3 in DIR])],
    [ax_lib_fftw3_with_fftw3="$withval"], [ax_lib_fftw3_with_fftw3=])

  AC_ARG_WITH(fftw3-prefix,
    [AC_HELP_STRING([--with-fftw3-prefix=PREFIX], [compile with fftw3 prefix PREFIX])],
    [ax_lib_fftw3_with_fftw3_prefix="$withval"], [ax_lib_fftw3_with_fftw3_prefix=])

  AC_ARG_WITH(fftw3-includedir,
    [AC_HELP_STRING([--with-fftw3-includedir=DIR], [compile with fftw3 include directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_fftw3_with_fftw3_inc_dir="$withval"], [ax_lib_fftw3_with_fftw3_inc_dir=])

  AC_ARG_WITH(fftw3-libdir,
    [AC_HELP_STRING([--with-fftw3-libdir=DIR], [compile with fftw3 library directory DIR (DIR can be a whitespace-separated list of directories)])],
    [ax_lib_fftw3_with_fftw3_lib_dir="$withval"], [ax_lib_fftw3_with_fftw3_lib_dir=])
])


AC_DEFUN([AX_LIB_FFTW3_CHECK],[
  _AX_LIB_FFTW3_CHECK([],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_FFTW3F_CHECK],[
  _AX_LIB_FFTW3_CHECK([f],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_FFTW3L_CHECK],[
  _AX_LIB_FFTW3_CHECK([l],[$1],[$2],[$3],[$4])
])

AC_DEFUN([AX_LIB_FFTW3Q_CHECK],[
  _AX_LIB_FFTW3_CHECK([q],[$1],[$2],[$3],[$4])
])


AC_DEFUN([_AX_LIB_FFTW3_CHECK],[

  ax_type_suffix="$1"
  ax_with_fftw3=m4_ifnblank($2,"$2","$ax_lib_fftw3_with_fftw3")
  ax_with_fftw3_prefix=m4_ifnblank($3,"$3","$ax_lib_fftw3_with_fftw3_prefix")
  ax_with_fftw3_inc_dir=m4_ifnblank($4,"$4","$ax_lib_fftw3_with_fftw3_inc_dir")
  ax_with_fftw3_lib_dir=m4_ifnblank($5,"$5","$ax_lib_fftw3_with_fftw3_lib_dir")

#  echo "with_fftw3: $ax_with_fftw3"
#  echo "with_fftw3_prefix: $ax_with_fftw3_prefix"
#  echo "with_fftw3_lib_dir: $ax_with_fftw3_lib_dir"
#  echo "with_fftw3_inc_dir: $ax_with_fftw3_inc_dir"

  fftw3_CPPFLAGS=
  fftw3_LDFLAGS=
  fftw3_PREFIX=
  fftw3_LIBS=
  fftw3_threads_LIBS=
  fftw3_mpi_LIBS=

  if test "x$ax_with_fftw3" != x ; then
    if test "x$ax_with_fftw3_inc_dir" = x ; then
      ax_with_fftw3_inc_dir="$ax_with_fftw3/include"
    fi
    if test "x$ax_with_fftw3_lib_dir" = x ; then 
      ax_with_fftw3_lib_dir="$ax_with_fftw3/lib"
    fi
  fi

  if test "x$ax_with_fftw3_inc_dir" != x ; then
    for ax_d in $ax_with_fftw3_inc_dir ; do
      if test "x$ax_with_fftw3" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The include directory '$ax_d' does not exist.])])
      fi
      fftw3_CPPFLAGS="$fftw3_CPPFLAGS -I$ax_d"
    done
    fftw3_CPPFLAGS="${fftw3_CPPFLAGS# }"
  fi

  if test "x$ax_with_fftw3_lib_dir" != x ; then 
    for ax_d in $ax_with_fftw3_lib_dir ; do
      if test "x$ax_with_fftw3" != xbuilt ; then
        AX_CHECK_DIR([$ax_d],[],[
          AC_MSG_WARN([The library directory '$ax_d' does not exist.])])
      fi
      fftw3_LDFLAGS="$fftw3_LDFLAGS -L$ax_d"
    done
    fftw3_LDFLAGS="${fftw3_LDFLAGS# }"
  fi

  AC_LANG_PUSH([C])

  saved_CFLAGS="$CFLAGS"
  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  CPPFLAGS="$CPPFLAGS $fftw3_CPPFLAGS"
  LDFLAGS="$LDFLAGS $fftw3_LDFLAGS"

  ax_lib_fftw3=no
  ax_lib_fftw3_threads=no
  ax_lib_fftw3_mpi=no

  fftw3_PREFIX="$ax_with_fftw3_prefix"

  if test "x$ax_with_fftw3" = xbuilt ; then
#    AC_MSG_NOTICE([fftw is built!])
    if test -z "$ax_with_fftw3_inc_dir" || test -z "$ax_with_fftw3_lib_dir"; then
      AC_MSG_ERROR([for an in-tree fftw3, set include and library directory with --with-fftw3-include-dir and --with-fftw3-lib-dir])
    fi
    ax_lib_fftw3=yes
    fftw3_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}"

    ax_lib_fftw3_threads=yes
    case $ac_configure_args in
      *--enable-threads*)
        fftw3_threads_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_threads"
        ;;
      *)
        ax_lib_fftw3_threads=no
        ;;
    esac
    
    ax_lib_fftw3_mpi=yes
    case $ac_configure_args in
      *--enable-mpi*)
        fftw3_mpi_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_mpi"
        ;;
      *)
        ax_lib_fftw3_mpi=no
        ;;
    esac
  else
    # Check if header is present and usable.
    AC_CHECK_HEADER([fftw3.h], [ax_lib_fftw3=yes; ax_lib_fftw3_threads=yes])
    AC_CHECK_HEADER([fftw3-mpi.h], [ax_lib_fftw3_mpi=yes])

    if test "x$ax_lib_fftw3" = xyes ; then
      saved_LIBS="$LIBS"

      AC_CHECK_LIB([${ax_with_fftw3_prefix}fftw3${ax_type_suffix}], [${ax_with_fftw3_prefix}fftw${ax_type_suffix}_execute], [], [ax_lib_fftw3=no])
      test "x${ax_lib_fftw3}" != xno && fftw3_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}"
      if test "x$ax_lib_fftw3_threads" = xyes ; then
        AC_CHECK_LIB([${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_threads], [${ax_with_fftw3_prefix}fftw${ax_type_suffix}_init_threads], [], [ax_lib_fftw3_threads=no])
        test "x${ax_lib_fftw3_threads}" != xno && fftw3_threads_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_threads"
      fi
      if test "x$ax_lib_fftw3_mpi" = xyes ; then
        AC_CHECK_LIB([${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_mpi], [${ax_with_fftw3_prefix}fftw${ax_type_suffix}_mpi_init], [], [ax_lib_fftw3_mpi=no])
        test "x${ax_lib_fftw3_mpi}" != xno && fftw3_mpi_LIBS="-l${ax_with_fftw3_prefix}fftw3${ax_type_suffix}_mpi"
      fi
    
      LIBS="$saved_LIBS"
    fi
  fi

  # Restore saved flags.
  CFLAGS="$saved_CFLAGS"
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"

  AC_LANG_POP([C])

  AC_SUBST(fftw3_CPPFLAGS)
  AC_SUBST(fftw3_LDFLAGS)
  AC_SUBST(fftw3_PREFIX)
  AC_SUBST(fftw3_LIBS)
  AC_SUBST(fftw3_threads_LIBS)
  AC_SUBST(fftw3_mpi_LIBS)

  if test "x$ax_lib_fftw3" = xyes ; then
    AC_DEFINE_UNQUOTED([FFTW_PREFIX],[$fftw3_PREFIX],[Define to the prefix of the namespace of the FFTW library.])
  fi
])
