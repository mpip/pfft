# _AX_FC_INTEGER_4_FLAGS
# ----------------------
# Flags to set INTEGER kind to 4 bytes.
AC_DEFUN([_AX_FC_INTEGER_4_FLAGS],
[-fdefault-integer-4 -qintsize=4 "-integer-size 32" -CcdII4 "-s integer32" -xtypemap=integer:32 -i4 +i4 -nlong])

# _AX_FC_INTEGER_8_FLAGS
# ----------------------
# Flags to set INTEGER kind to 8 bytes.
AC_DEFUN([_AX_FC_INTEGER_8_FLAGS],
[-fdefault-integer-8 -qintsize=8 "-integer-size 64" -CcdII8 "-s integer64" -xtypemap=integer:64 -i8 +i8 -long])

# AX_FC_INTEGER_KIND(KIND, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------------
AC_DEFUN([AX_FC_INTEGER_KIND],
[m4_case([$1], [4], [], [8], [],
   [m4_fatal([$0: KIND argument must be 4 or 8])])dnl
AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK(
[for Fortran flag needed to enable $1 byte default integer kinds],
	       [ax_cv_fc_integer_kind_$1],
[ax_cv_fc_integer_kind_$1=unknown
ax_fc_integer_kind_FCFLAGS_save=$FCFLAGS
for ax_flag in none _AX_FC_INTEGER_]$1[_FLAGS
do
  test "x$ax_flag" != xnone && FCFLAGS="$ax_fc_integer_kind_FCFLAGS_save $ax_flag"
  AX_CHECK_KINDOF([integer], [silent])
  if test "$ax_check_kindof" = $1; then
    ax_cv_fc_integer_kind_$1=$ax_flag
    break
  fi
done
FCFLAGS=$ax_fc_integer_kind_FCFLAGS_save
])
if test "x$ax_cv_fc_integer_kind_$1" = xunknown; then
  m4_default([$3], [:])
else
  m4_default([$2], [:])
fi
AC_LANG_POP([Fortran])dnl
])
