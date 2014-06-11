# Copyright (C) 2011 The ScaFaCoS project
#  
# This file is part of ScaFaCoS.
#  
# ScaFaCoS is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
#  ScaFaCoS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#


# add configure option --enable-dist
AC_DEFUN_ONCE([AX_FCS_ENABLE_DIST_ARG],[
AC_ARG_ENABLE(dist,
  AC_HELP_STRING([--enable-dist],[Enable everything to make a complete tarball with 'make dist'.]),,[enable_dist=no])
])


# add configure option --enable-single-lib
AC_DEFUN_ONCE([AX_FCS_ENABLE_SINGLE_LIB_ARG],[
AC_ARG_ENABLE([single-lib],
  AC_HELP_STRING([--enable-single-lib],[Create only a single libfcs.a file containing the whole library.]),,[enable_single_lib=no])
])

AC_DEFUN_ONCE([AX_FCS_ENABLE_SINGLE_LIB],[
AM_CONDITIONAL([ENABLE_SINGLE_LIB],[test "x$enable_single_lib" = xyes])
])


# add configure option --disable-library-install
AC_DEFUN_ONCE([AX_FCS_DISABLE_LIBRARY_INSTALL_ARG],[
m4_ifdef([AX_DISABLE],[AX_DISABLE([library-install])])
])

AC_DEFUN_ONCE([AX_FCS_DISABLE_LIBRARY_INSTALL],[
AM_CONDITIONAL([ENABLE_LIBRARY_INSTALL],[test "x$ax_disable_library_install" != xyes])
])


# List of all switches and modules.
m4_define([fcs_switches], [info debug timing])
m4_define([fcs_modules], [near gridsort resort sl_fmm pnfft])

# add configure options --enable-fcs-info, --enable-fcs-debug, and --enable-fcs-timing
AC_DEFUN_ONCE([AX_FCS_ENABLE_INFODEBUG_ARG],[
AC_FOREACH([switch], fcs_switches,[
  AC_MSG_CHECKING([whether to enable ]switch) 
  AC_ARG_ENABLE(fcs-switch,
    AC_HELP_STRING([--enable-fcs-]switch,[enable ]switch[ output]),,[enable_fcs_]switch[=no])
  AC_MSG_RESULT([$enable_fcs_]switch)
  ])
])

AC_DEFUN_ONCE([AX_FCS_ENABLE_INFODEBUG],[
# Let debug enable info automatically.
test "x$enable_fcs_debug" != xno && enable_fcs_info="$enable_fcs_debug"
AC_FOREACH([switch], fcs_switches,[
if test "x$enable_fcs_[]switch[]" != xno ; then
  AC_DEFINE([FCS_ENABLE_]m4_toupper(switch), [1], [Define if ]switch[ output is enabled.])
  AC_FOREACH([module], fcs_modules,[
    case $enable_fcs_]switch[ in
      *module*)
        AC_DEFINE([FCS_ENABLE_]m4_toupper(switch)[_]m4_toupper(module), [1], [Define if ]switch[ output of ]module[ is enabled.])
        ;;
    esac
  ])
fi
])
])


# add configure option --disable-fcs-fortran
AC_DEFUN_ONCE([AX_FCS_DISABLE_FORTRAN_ARG],[
AC_ARG_ENABLE([fcs-fortran],
  AC_HELP_STRING([--disable-fcs-fortran],[Disable Fortran interface and tests.]),,[enable_fcs_fortran=yes])
])

AC_DEFUN_ONCE([AX_FCS_DISABLE_FORTRAN],[
if test "x$enable_fcs_fortran" = xyes ; then
  AC_DEFINE([FCS_ENABLE_FORTRAN], [1],[Define if Fortran interface is enabled.])
  use_fcs_fortran=yes
else
  use_fcs_fortran=no
fi
AM_CONDITIONAL([ENABLE_FORTRAN_INTERFACE],[test "x$use_fcs_fortran" = xyes])
])


# add configure option --enable-fcs-extended-tests
AC_DEFUN_ONCE([AX_FCS_ENABLE_EXTENDED_TESTS_ARG],[
AC_ARG_ENABLE([fcs-extended-tests],
  AC_HELP_STRING([--enable-fcs-extended-tests],[Enable extended tests of the FCS library and its solvers and subpackages.]),,[enable_fcs_extended_tests=no])
])


# AX_FCS_CONFIG_SOLVER_ARGS
# -------------------------
# Set up FCS library stuff arguments in solver configure.
AC_DEFUN_ONCE([AX_FCS_CONFIG_SOLVER_ARGS],[
AX_FCS_ENABLE_DIST_ARG
AX_FCS_DISABLE_LIBRARY_INSTALL_ARG
AX_FCS_ENABLE_INFODEBUG_ARG
AX_FCS_TYPES_ARGS
])


# AX_FCS_CONFIG_SOLVER
# --------------------
# Set up FCS library stuff in solver configure.
AC_DEFUN([AX_FCS_CONFIG_SOLVER],[
AX_FCS_CONFIG_SOLVER_ARGS
AX_FCS_DISABLE_LIBRARY_INSTALL
AX_FCS_ENABLE_INFODEBUG
AX_FCS_TYPES
])


# AX_FCS_CONFIG_TOP_ARGS
# ----------------------
# Set up FCS library stuff arguments in top-level configure.
AC_DEFUN_ONCE([AX_FCS_CONFIG_TOP_ARGS],[
AX_FCS_ENABLE_DIST_ARG
AX_FCS_ENABLE_SINGLE_LIB_ARG
AX_FCS_ENABLE_INFODEBUG_ARG
AX_FCS_DISABLE_FORTRAN_ARG
AX_FCS_ENABLE_EXTENDED_TESTS_ARG
AX_FCS_TYPES_ARGS
])


# AX_FCS_CONFIG_TOP
# -----------------
# Set up FCS library stuff in top-level configure.
AC_DEFUN([AX_FCS_CONFIG_TOP],[
AX_FCS_CONFIG_TOP_ARGS
AX_FCS_ENABLE_SINGLE_LIB
AX_FCS_ENABLE_INFODEBUG
AX_FCS_DISABLE_FORTRAN
AX_FCS_TYPES
])
