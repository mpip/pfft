# This macro has been submitted for inclusion into the Autoconf 2.69 tree.
# Use the upstream version when it is defined.
m4_version_prereq([2.69], [],
[

# AC_FC_PP_DEFINE([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -------------------------------------------------------------------
# Find a flag to specify defines for preprocessed Fortran.  Not all
# Fortran compilers use -D.  Substitutes FC_DEFINE with the result and
# calls ACTION-IF-SUCCESS (defaults to nothing) if successful, and
# ACTION-IF-FAILURE (defaults to failing with an error message) if not.
#
# Known flags:
# IBM: -WF,-D
# Lahey/Fujitsu: -Wp,-D     older versions???
# f2c: -D or -Wc,-D
# others: -D
AC_DEFUN([AC_FC_PP_DEFINE],
[AC_LANG_PUSH([Fortran])dnl
ac_fc_pp_define_srcext_save=$ac_fc_srcext
AC_FC_SRCEXT([F])
AC_CACHE_CHECK([how to define symbols for preprocessed Fortran],
  [ac_cv_fc_pp_define],
[ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_cv_fc_pp_define=unknown
ac_fc_pp_define_FCFLAGS_save=$FCFLAGS
for ac_flag in -D -WF,-D -Wp,-D -Wc,-D
do
  FCFLAGS="$ac_fc_pp_define_FCFLAGS_save ${ac_flag}FOOBAR ${ac_flag}ZORK=42"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#ifndef FOOBAR
      choke me
#endif
#if ZORK != 42
      choke me
#endif]])],
    [ac_cv_fc_pp_define=$ac_flag; break])
done
FCFLAGS=$ac_fc_pp_define_FCFLAGS_save
])
ac_fc_srcext=$ac_fc_pp_define_srcext_save
if test "x$ac_cv_fc_pp_define" = xunknown; then
  FC_DEFINE=
  m4_default([$2],
	     [AC_MSG_ERROR([Fortran does not allow to define preprocessor symbols], 77)])
else
  FC_DEFINE=$ac_cv_fc_pp_define
  $1
fi
AC_SUBST([FC_DEFINE])dnl
AC_LANG_POP([Fortran])dnl
])

])dnl m4_version_prereq([2.69])
