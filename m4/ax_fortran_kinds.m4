# AX_FORTRAN_KINDS
# ----------------

# AX_CHECK_KINDOF(FTYPE[,VOLATILE])
# ---------------------------------
# Determine the kind of Fortran type FTYPE and assign the value to
# ax_check_kindof. If FTYPE is not a shell variable (i.e., contains no $) and
# VOLATILE is empty, then the same behavior as AC_CHECK_SIZEOF is used:
# KINDOF_<type> is defined to be the kind of FTYPE and the result is cached in
# ax_cv_kindof_<type>. Otherwise, no C preprocessor macro is defined and the
# result is not cached. If VOLATILE is "silent", then no messages with
# AC_MSG_... are shown.

AC_DEFUN([AX_CHECK_KINDOF],[
ax_check_kindof=
m4_if($1$2,m4_translit([[$1]],[$],[_]),
 [m4_define([sub_id],m4_translit(m4_strip([[$1]]),[()=, ],[_____]))
  m4_define([var_id],[ax_cv_kindof_]AS_TR_SH(sub_id))
  m4_define([cpp_id],KINDOF_[]AS_TR_CPP(sub_id))
#  echo "sub: sub_id"
#  echo "var: var_id"
#  echo "cpp: cpp_id"
  AC_CACHE_CHECK([kind of $1],[var_id],
   [AX_FORTRAN_KINDS_CHECK_RUN([ax_check_kindof],$1)
    if test "${ax_check_kindof}" = cross -o "${ax_check_kindof}" = failed ; then
      AX_FORTRAN_KINDS_CHECK_COMPILE([ax_check_kindof],$1)
    fi
    var_id="${ax_check_kindof}"
  ])
  ax_check_kindof="${var_id}"
  if test "x${ax_check_kindof}" != x ; then
    AC_DEFINE_UNQUOTED([cpp_id],${ax_check_kindof},[The kind of '$1'.])
  fi
 ],
 [m4_if($2,[silent],,AC_MSG_CHECKING([kind of $1 (volatile)]))
  AX_FORTRAN_KINDS_CHECK_RUN([ax_check_kindof],$1)
  if test "${ax_check_kindof}" = cross -o "${ax_check_kindof}" = failed ; then
    AX_FORTRAN_KINDS_CHECK_COMPILE([ax_check_kindof],$1)
  fi
  m4_if($2,[silent],,AC_MSG_RESULT([$ax_check_kindof]))
])
])


# AX_FORTRAN_KINDS_CHECK_RUN(VAR, FTYPE)
# --------------------------------------
# Determine the kind of Fortran data type FTYPE (by running Fortran code) and
# assign the value to VAR. Set VAR to "failed" or "cross" if the code fails or
# cross compilation is enabled.

AC_DEFUN([AX_FORTRAN_KINDS_CHECK_RUN],[
AC_LANG_PUSH([Fortran])dnl
AC_RUN_IFELSE(
 [AC_LANG_PROGRAM([],[[
      $2 :: v
      open(90,file='conftest.val')
      write (90,'(i2)') kind(v)
      close(90)
]])],
 [read $1 < conftest.val],[$1=failed],[$1=cross])dnl
AC_LANG_POP([Fortran])dnl
#echo "run: $ax_check_kindof"
])


# AX_FORTRAN_KINDS_CHECK_COMPILE(VAR, FTYPE[, KINDS])
# ---------------------------------------------------
# Determine the kind of Fortran data type FTYPE (by compiling Fortran code) and
# assign the value to VAR. Set VAR to "failed" if the kind could not be
# determined. KINDS can be a list of kind values to check. If KINDS is empty,
# then kind values "2 4 8 16" are used. The test requires that the Fortran
# compiler fails if arrays of different size are assigned to each other
# (e.g., "a(1:2) = b(1:3)" should fail!).

AC_DEFUN([AX_FORTRAN_KINDS_CHECK_COMPILE],[
ax_kind=
if test "x$3" != x ; then
  ax_kind_list="$3"
else
  ax_kind_list="2 4 8 16"
fi
AC_LANG_PUSH([Fortran])dnl
for ax_k in ${ax_kind_list} ; do
  AC_COMPILE_IFELSE(
   [AC_LANG_PROGRAM([],[[
      implicit none
      $2 :: dt
      integer, parameter :: dk = kind(dt)
      integer :: dka(dk), tka(${ax_k})
      dka(1:dk) = tka(1:${ax_k})
]])],
   [ax_kind="${ax_k}"])
  test "x${ax_kind}" != x && break
done
AC_LANG_POP([Fortran])dnl
if test "x${ax_kind}" != x ; then
  $1="${ax_kind}"
else
  $1=failed
fi
#echo "compile: $ax_check_kindof"
])


# AX_FORTRAN_KINDS_DEFAULT
# ------------------------
# Determine the default kinds of integer and real types in Fortran. Assign the
# results to ax_default_fortran_integer_kind and ax_default_fortran_real_kind
# and define DEFAULT_FORTRAN_INTEGER_KIND and DEFAULT_FORTRAN_REAL_KIND to the
# default kinds. Both results are cached in ax_cv_default_fortran_kinds.

AC_DEFUN([AX_FORTRAN_KINDS_DEFAULT],[
AC_CACHE_CHECK(
 [default integer / real kind],
 [ax_cv_default_fortran_kinds],
 [if test "x${ax_cv_kindof_integer}" != x ; then
    ax_cv_default_fortran_kinds="${ax_cv_kindof_integer}"
  else
    AX_CHECK_KINDOF([integer],[silent])
    ax_cv_default_fortran_kinds="${ax_check_kindof}"
  fi
  if test "x${ax_cv_kindof_real}" != x ; then
    ax_cv_default_fortran_kinds="${ax_cv_default_fortran_kinds} / ${ax_cv_kindof_real}"
  else
    AX_CHECK_KINDOF([real],[silent])
    ax_cv_default_fortran_kinds="${ax_cv_default_fortran_kinds} / ${ax_check_kindof}"
  fi
])
ax_default_fortran_integer_kind="${ax_cv_default_fortran_kinds%% /*}"
AC_DEFINE_UNQUOTED([DEFAULT_FORTRAN_INTEGER_KIND],[${ax_default_fortran_integer_kind}],[Default kind of Fortran integer data type.])
ax_default_fortran_real_kind="${ax_cv_default_fortran_kinds##*/ }"
AC_DEFINE_UNQUOTED([DEFAULT_FORTRAN_REAL_KIND],[${ax_default_fortran_real_kind}],[Default kind of Fortran real data type.])
])


# AX_FORTRAN_C2F_TYPE(VAR, CTYPE, FTYPES)
# ---------------------------------------
# Check which Fortran type in FTYPES (comma-separated list) matches C type CTYPE
# and assign the value to VAR.

AC_DEFUN([AX_FORTRAN_C2F_TYPE],[
ax_c=failed
m4_foreach(ct,m4_split([char,short,int,long,long long,float,double,long double],[,]),
 [if test "$2" = "ct" ; then
    AC_CHECK_SIZEOF(m4_strip(ct))
    ax_c="${ac_cv_sizeof_[]AS_TR_SH(m4_strip(ct))}"
  fi
 ])
ax_ft=failed
if test "x${ax_c}" != xfailed ; then
  m4_foreach(ft,m4_split([$3],[,]),
   [if test "${ax_ft}" = failed ; then
      AX_CHECK_KINDOF(ft)
      test "${ax_c}" = "${ax_check_kindof}" && ax_ft="m4_strip(ft)"
    fi
   ])
fi
$1="${ax_ft}"
])


# AX_FORTRAN_F2C_TYPE(VAR, FTYPE, CTYPES)
# ---------------------------------------
# Check which C type in CTYPES (comma-separated list) matches Fortran type FTYPE
# and assign the value to VAR.

AC_DEFUN([AX_FORTRAN_F2C_TYPE],[
AX_CHECK_KINDOF($2)
ax_f="${ax_check_kindof}"
ax_ct=failed
if test "x${ax_f}" != xfailed ; then
  m4_foreach(ct,m4_split([$3],[,]),
   [if test "${ax_ct}" = failed ; then
      AC_CHECK_SIZEOF(m4_strip(ct))
      test "${ax_f}" = "${ac_cv_sizeof_[]AS_TR_SH(m4_strip(ct))}" && ax_ct="m4_strip(ct)"
    fi
   ])
fi
$1="${ax_ct}"
])


# AX_FORTRAN_C2F_KIND_ISOC(VAR, CTYPE)
# ------------------------------------
# Determine the ISO C kind of C type CTYPE and assign the value to VAR.

AC_DEFUN([AX_FORTRAN_C2F_KIND_ISOC],[
  case $2 in
    "short") $1="c_short" ;;
    "int") $1="c_int" ;;
    "long") $1="c_long" ;;
    "long long") $1="c_long_long" ;;
    "float") $1="c_float" ;;
    "double") $1="c_double" ;;
    "long double") $1="c_long_double" ;;
    *) $1=failed ;;
  esac
])
