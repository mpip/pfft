# AX_FC_NO_RANGE_CHECK([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# --------------------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) disable range checking on
# results of simplification of constant expressions during compilation, and
# add it to FCFLAGS.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful and
# ACTION-IF-FAILURE (defaults to nothing) if not.
#
# The known flags are:
# gfortran: -fno-range-check
# HP: +check=truncate:none

AC_DEFUN([AX_FC_NO_RANGE_CHECK],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK(
[for Fortran flag needed to turn off range checking],
	       [ax_cv_fc_no_range_check],
[ax_cv_fc_no_range_check=unknown
ax_fc_no_range_check_FCFLAGS_save=$FCFLAGS
for ax_flag in none -fno-range-check +check=truncate:none
do
  test "x$ax_flag" != xnone && FCFLAGS="$ax_fc_no_range_check_FCFLAGS_save $ax_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
      INTEGER (kind=4) :: i
      i = -2147483648]])],
		    [ax_cv_fc_no_range_check=$ax_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ax_fc_no_range_check_FCFLAGS_save
])
if test "x$ax_cv_fc_no_range_check" = xunknown; then
  m4_default([$3], [:])
else
  if test "x$ax_cv_fc_no_range_check" != xnone; then
    FCFLAGS="$FCFLAGS $ax_cv_fc_no_range_check"
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# AX_FC_NO_RANGE_CHECK
