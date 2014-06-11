# This macro has been submitted for inclusion into the Autoconf 2.69 tree.
# Use the upstream version when it is defined.
m4_version_prereq([2.69], [],
[

# AC_FC_PP_SRCEXT(EXT, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# --------------------------------------------------------------
# Like AC_FC_SRCEXT, set the source-code extension used in Fortran (FC) tests
# to EXT (which defaults to f).  Also, look for any necessary additional
# FCFLAGS needed to allow this extension for preprocessed Fortran, and store
# them in the output variable FCFLAGS_<EXT> (e.g. FCFLAGS_f90 for EXT=f90).
# If successful, call ACTION-IF-SUCCESS.  If unable to compile preprocessed
# source code with EXT, call ACTION-IF-FAILURE, which defaults to failing with
# an error message.
#
# Some compilers allow preprocessing with either a Fortran preprocessor or
# with the C preprocessor (cpp).  Prefer the Fortran preprocessor, to deal
# correctly with continuation lines, `//' (not a comment), and preserve white
# space (for fixed form).
#
# (The flags for the current source-code extension, if any, are stored in
# $ac_fcflags_srcext and used automatically in subsequent autoconf tests.)
#
# For ordinary extensions like f90, etcetera, the modified FCFLAGS
# are needed for IBM's xlf*.  Also, for Intel's ifort compiler, the
# $FCFLAGS_<EXT> variable *must* go immediately before the source file on the
# command line, unlike other $FCFLAGS.  Ugh.
#
# gfortran 4.4.3 does not fail upon nonexistent header file, nor upon
# the standard Autoconf fail string `choke me'.  Go figure.
#
# Known extensions that enable preprocessing by default, and flags to force it:
# GNU: .F .F90 .F95 .F03 .F08, -cpp for most others,
#      -x f77-cpp-input for .f77 .F77; -x f95-cpp-input for gfortran < 4.4
# SGI: .F .F90, -ftpp or -cpp for .f .f90
# SUN: .F .F95, -fpp for others; -xpp={fpp,cpp} for preprocessor selection
# IBM: .F .F77 .F90 .F95 .F03, -qsuffix=cpp=EXT for extension .EXT to invoke cpp
#      -WF,-qnofpp -WF,-qfpp=comment:linecont:nocomment:nolinecont
#      -WF,-qlanglvl=classic or not -qnoescape (trigraph problems)
# HP:  .F, +cpp={yes|no|default} use cpp, -cpp
# PGI: -Mpreprocess
# Absoft: .F .FOR .F90 .F95, -cpp for others
# Cray: .F .F90 .FTN, -e Z for others
# Intel: .F .F90, -fpp for others, but except for .f and .f90, -Tf may also be
#        needed right before the source file name
# PathScale: .F .F90 .F95, -ftpp or -cpp for .f .f90 .f95
# Lahey: .F .FOR .F90 .F95, -Cpp
AC_DEFUN([AC_FC_PP_SRCEXT],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile preprocessed .$1 files],
		ac_cv_fc_pp_srcext_$1,
[ac_ext=$1
# Also adjust FCFLAGS_SRCEXT, for compatibility with Autoconf 2.59.
ac_FCFLAGS_SRCEXT_save=$FCFLAGS_SRCEXT
ac_fcflags_pp_srcext_save=$ac_fcflags_srcext
ac_fcflags_srcext=
FCFLAGS_SRCEXT=
ac_cv_fc_pp_srcext_$1=unknown
case $ac_ext in #(
  [[fF]]77) ac_try=f77-cpp-input;; #(
  *) ac_try=f95-cpp-input;;
esac
for ac_flag in none -ftpp -fpp -Tf "-fpp -Tf" -xpp=fpp -Mpreprocess "-e Z" \
               -cpp -xpp=cpp -qsuffix=cpp=$1 "-x $ac_try" +cpp -Cpp; do
  if test "x$ac_flag" != xnone; then
    ac_fcflags_srcext=$ac_flag
    FCFLAGS_SRCEXT=$ac_flag
  fi
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 0
#include <ac_nonexistent.h>
this is not correct fortran
#endif]])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 1
#include <ac_nonexistent.h>
this is not correct fortran
#endif]])],
       [],
       [ac_cv_fc_pp_srcext_$1=$ac_flag; break])])
done
rm -f conftest.$ac_objext conftest.$1
ac_fcflags_srcext=$ac_fcflags_pp_srcext_save
FCFLAGS_SRCEXT=$ac_FCFLAGS_SRCEXT_save=
])
if test "x$ac_cv_fc_pp_srcext_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran could not compile preprocessed .$1 files])])
else
  ac_fc_srcext=$1
  if test "x$ac_cv_fc_pp_srcext_$1" = xnone; then
    ac_fcflags_srcext=""
    FCFLAGS_SRCEXT=
    FCFLAGS_[]$1[]=""
  else
    ac_fcflags_srcext=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_SRCEXT=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_[]$1[]=$ac_cv_fc_pp_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
AC_LANG_POP(Fortran)dnl
])# AC_FC_PP_SRCEXT

])dnl m4_version_prereq([2.69])
