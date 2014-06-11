# This is a modified version of AC_PROG_CXXCPP from Autoconf
# changed to cater to Fortran needs.  Some additional helper
# macros are included below as well.

# ACX_PROG_FCCPP
# --------------
# Find a preprocessor suitable for preprocessing Fortran.
# Note that the file name extension may be important for the
# preprocessor, because // may be treated as comment in non-Fortran
# mode of `$FC -E'.
AC_DEFUN([ACX_PROG_FCCPP],
[AC_REQUIRE([AC_PROG_CPP])dnl
AC_ARG_VAR([FCCPP], [Preprocessor for Fortran])dnl
AC_ARG_VAR([FCCPPFLAGS],
	   [Preprocessor flags for the Fortran preprocessor])
AC_LANG_PUSH([Fortran])dnl
AC_MSG_CHECKING([how to preprocess Fortran])
if test -z "$FCCPP"; then
  AC_CACHE_VAL([ac_cv_prog_FCCPP],
    [# gfortran needs -cpp.
     # gfortran -cpp -E keeps // comments if passed a .f or .f90 file.
     # Don't try -traditional-cpp, breaks var args
     for FCCPP in "cpp -C -P" cpp \
                  "/lib/cpp -C -P" /lib/cpp \
		  "$FC -cpp -E" "$FC -E" "$CPP -C -P" "$CPP" "$CPP -x c -C -P"
     do
       _AC_PROG_FCPREPROC_WORKS_IFELSE([ac_cv_prog_FCCPP=$FCCPP; break])
     done
    ])
  FCCPP=$ac_cv_prog_FCCPP
fi
AC_MSG_RESULT([$FCCPP])
AC_SUBST([FCCPP])dnl
# Default FCCPPFLAGS to CPPFLAGS.
: ${FCCPPFLAGS=$CPPFLAGS}
if test -z "$FCCPP"; then
   AC_MSG_FAILURE([
  No viable Fortran preprocessor found, please set FCCPP environment variable. The
  Fortran preprocessor has to support variadic macros and put these on a single line.
  A gcc-compatible C preprocessor works, the clang preprocessor doesn't.])
fi
_AC_PROG_FCPREPROC_WORKS_IFELSE([],
          [AC_MSG_FAILURE([Fortran preprocessor "$FCCPP" fails sanity check. For a gcc-like cpp,
you may have to add flags "-C -P".])])
AC_LANG_POP([Fortran])dnl
])

# _AC_PROG_FCPREPROC_WORKS_IFELSE([IF-WORKS], [IF-NOT])
# -----------------------------------------------------
# Helper macro to find out whether Fortran preprocessing works ok.
AC_DEFUN([_AC_PROG_FCPREPROC_WORKS_IFELSE],
[
# for checking the variadic macro output
AC_REQUIRE([AC_PROG_EGREP])
ac_cpp_SAVE=$ac_cpp
ac_cpp='$FCCPP $FCCPPFLAGS'
: >conftest.h
for ac_[]_AC_LANG_ABBREV[]_preproc_warn_flag in '' yes
do
  ac_preproc_fail=false
  ac_preproc_pass=false
  # For Fortran, "//" may not be discarded because it is part of the
  # valid syntax.
  _AC_FCPREPROC_IFELSE([AC_LANG_SOURCE([[
@%:@ include "conftest.h"
		     keep//this]])],
                     [fc_grep='grep keep//this conftest.i >/dev/null'
		      if AC_TRY_EVAL([fc_grep])
		      then
			ac_preproc_pass=:
		      else
			continue
		      fi],
                     [# Broken: fails on valid input.
continue])

  # pepc and fmm require variadic macros to be crammed into a single line,
  # regardless of how they were formatted. gcc does this for now, clang
  # doesn't. Better check...
  _AC_FCPREPROC_IFELSE([AC_LANG_SOURCE([[
@%:@ include "conftest.h"
/* intentional indentation, clang also doesn't like spaces before define */
 #define varargs(...) __VA_ARGS__
 varargs(this_is,
         on_one_line)]])],
                     [fc_egrep='$EGREP this_is.*on_one_line conftest.i >/dev/null'
		      if AC_TRY_EVAL([fc_egrep])
		      then
			ac_preproc_pass=:
		      else
			continue
		      fi],
                     [# Broken: fails on valid input.
continue])

  # OK, works on sane cases.  Now check whether nonexistent headers
  # can be detected and how.
  _AC_FCPREPROC_IFELSE([AC_LANG_SOURCE([[@%:@include <ac_nonexistent.h>]])],
                     [# Broken: success on invalid input.
continue],
                     [# Passes the second test.
ac_preproc_fail=:
break])

done
ac_cpp=$ac_cpp_SAVE
# Because of `break', _AC_PREPROC_IFELSE's cleaning code was skipped.
rm -f conftest.h conftest.err conftest.$ac_ext
AS_IF([$ac_preproc_pass && $ac_preproc_fail], [$1], [$2])
])# _AC_PROG_FCPREPROC_WORKS_IFELSE


# _AC_FCPREPROC_IFELSE(PROGRAM, [ACTION-IF-TRUE], [ACTION-IF-FALSE])
# ----------------------------------------------------------------
# Try to preprocess PROGRAM.
#
# This macro can be used during the selection of a preprocessor.
# Run cpp and set ac_cpp_err to "yes" for an error, to
# "$ac_(c,cxx)_preproc_warn_flag" if there are warnings or to "" if
# neither warnings nor errors have been detected.  eval is necessary
# to expand ac_cpp.
AC_DEFUN([_AC_FCPREPROC_IFELSE],
[m4_ifvaln([$1], [AC_LANG_CONFTEST([$1])])dnl
if _AC_EVAL_STDERR([$ac_cpp conftest.$ac_ext]) > conftest.i; then
  if test -s conftest.err; then
    ac_cpp_err=$ac_[]_AC_LANG_ABBREV[]_preproc_warn_flag
    ac_cpp_err=$ac_cpp_err$ac_[]_AC_LANG_ABBREV[]_werror_flag
  else
    ac_cpp_err=
  fi
else
  ac_cpp_err=yes
fi
if test -z "$ac_cpp_err"; then
  m4_default([$2], :)
else
  _AC_MSG_LOG_CONFTEST
  $3
fi
rm -f conftest.err conftest.i m4_ifval([$1], [conftest.$ac_ext])[]dnl
])# _AC_FCPREPROC_IFELSE


# AC_FCPREPROC_IFELSE(PROGRAM, [ACTION-IF-TRUE], [ACTION-IF-FALSE])
# -----------------------------------------------------------------
# Try to preprocess PROGRAM.  Requires that the preprocessor for the
# current language was checked for, hence do not use this macro in macros
# looking for a preprocessor.
AC_DEFUN([AC_FCPREPROC_IFELSE],
[dnl AC_LANG_FCREPROC_REQUIRE()dnl
_AC_FCPREPROC_IFELSE($@)])

