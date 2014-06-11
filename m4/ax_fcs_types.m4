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
AC_DEFUN([AX_FCS_TYPE_INT_SET],[
AC_DEFINE([fcs_int], [$1], [Define to the integer type to use for FCS.])
AC_DEFINE([FCS_MPI_INT], [$2], [Define to the MPI datatype that corresponds to the integer type to use for FCS.])
AC_DEFINE([FCS_LMOD_INT], [$3], [Define to the printf length modifier that corresponds to the integer type to use for FCS.])
AC_DEFINE([FCS_CONV_INT], [$4], [Define to the scanf conversion that corresponds to the integer type to use for FCS.])
])


AC_DEFUN([AX_FCS_TYPE_INT],[
if test "$1" = "short" ; then
  AX_FCS_TYPE_INT_SET([short],[MPI_SHORT],["h"],["h"])
  AC_DEFINE([FCS_INT_IS_SHORT], [1], [Define whether fcs_int is short.])
elif test "$1" = "int" ; then
  AX_FCS_TYPE_INT_SET([int],[MPI_INT],[""])
  AC_DEFINE([FCS_INT_IS_INT], [1], [Define whether fcs_int is int.])
elif test "$1" = "long" ; then
  AX_FCS_TYPE_INT_SET([long],[MPI_LONG],["l"],["l"])
  AC_DEFINE([FCS_INT_IS_LONG], [1], [Define whether fcs_int is long.])
elif test "$1" = "long long" ; then
  AX_FCS_TYPE_INT_SET([long long],[MPI_LONG_LONG],["ll"],["ll"])
  AC_DEFINE([FCS_INT_IS_LONG_LONG], [1], [Define whether fcs_int is long long.])
else
  AC_MSG_ERROR([Datatype $1 is not a supported C integer type to use for FCS.])
fi
])


AC_DEFUN([AX_FCS_TYPE_FLOAT_SET],[
AC_DEFINE([fcs_float], [$1], [Define to the floating type to use for FCS.])
AC_DEFINE([FCS_MPI_FLOAT], [$2], [Define to the MPI datatype that corresponds to the floating type to use for FCS.])
AC_DEFINE([FCS_LMOD_FLOAT], [$3], [Define to the printf length modifier that corresponds to the floating type to use for FCS.])
AC_DEFINE([FCS_CONV_FLOAT], [$4], [Define to the scanf conversion that corresponds to the floating type to use for FCS.])
])


AC_DEFUN([AX_FCS_TYPE_FLOAT],[
if test "$1" = "float" ; then
  AX_FCS_TYPE_FLOAT_SET([float],[MPI_FLOAT],[""],[""])
  AC_DEFINE([FCS_FLOAT_IS_FLOAT], [1], [Define whether fcs_float is float.])
elif test "$1" = "double" ; then
  AX_FCS_TYPE_FLOAT_SET([double],[MPI_DOUBLE],[""],["l"])
  AC_DEFINE([FCS_FLOAT_IS_DOUBLE], [1], [Define whether fcs_float is double.])
elif test "$1" = "long double" ; then
  AX_FCS_TYPE_FLOAT_SET([long double],[MPI_LONG_DOUBLE],["L"],["L"])
  AC_DEFINE([FCS_FLOAT_IS_LONG_DOUBLE], [1], [Define whether fcs_float is long double.])
else
  AC_MSG_ERROR([Datatype $1 is not a supported C floating type to use for FCS.])
fi
])


AC_DEFUN([AX_FCS_TYPE_INTEGER_SET],[
AC_DEFINE([fcs_integer], [$1], [Define to the Fortran integer type to use for FCS, corresponding to fcs_int.])
AC_DEFINE([FCS_MPI_INTEGER], [$2], [Define to the MPI datatype that corresponds to the Fortran integer type to use for FCS.])
AC_DEFINE_UNQUOTED([fcs_integer_kind], [$3], [Define to the kind of fcs_integer.])
AC_DEFINE_UNQUOTED([fcs_integer_kind_isoc], [$4], [Define to the ISO C kind of fcs_integer.])
])


AC_DEFUN([AX_FCS_TYPE_INTEGER],[
AX_FORTRAN_C2F_KIND_ISOC([ax_kind_isoc],[$2])
if test "$1" = "integer" ; then
  AX_CHECK_KINDOF([integer])
  AX_FCS_TYPE_INTEGER_SET([integer],[MPI_INTEGER],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "integer*2" ; then
  AX_CHECK_KINDOF([integer*2])
  AX_FCS_TYPE_INTEGER_SET([integer*2],[MPI_INTEGER2],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "integer*4" ; then
  AX_CHECK_KINDOF([integer*4])
  AX_FCS_TYPE_INTEGER_SET([integer*4],[MPI_INTEGER4],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "integer*8" ; then
  AX_CHECK_KINDOF([integer*8])
  AX_FCS_TYPE_INTEGER_SET([integer*8],[MPI_INTEGER8],[${ax_check_kindof}],[${ax_kind_isoc}])
else
  AC_MSG_ERROR([Datatype $1 is not a supported Fortran integer type to use for FCS.])
fi
])


AC_DEFUN([AX_FCS_TYPE_REAL_SET],[
AC_DEFINE([fcs_real], [$1], [Define to the Fortran floating type to use for FCS, corresponding to fcs_float.])
AC_DEFINE([FCS_MPI_REAL], [$2], [Define to the MPI datatype that corresponds to the Fortran floating type to use for FCS.])
AC_DEFINE_UNQUOTED([fcs_real_kind], [$3], [Define to the kind of fcs_real.])
AC_DEFINE_UNQUOTED([fcs_real_kind_isoc], [$4], [Define to the ISO C kind of fcs_real.])
])


AC_DEFUN([AX_FCS_TYPE_REAL],[
AX_FORTRAN_C2F_KIND_ISOC([ax_kind_isoc],[$2])
if test "$1" = "real" ; then
  AX_CHECK_KINDOF([real])
  AX_FCS_TYPE_REAL_SET([real],[MPI_REAL],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "double precision" ; then
  AX_CHECK_KINDOF([double precision])
  AX_FCS_TYPE_REAL_SET([double precision],[MPI_DOUBLE_PRECISION],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "real*4" ; then
  AX_CHECK_KINDOF([real*4])
  AX_FCS_TYPE_REAL_SET([real*4],[MPI_REAL4],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "real*8" ; then
  AX_CHECK_KINDOF([real*8])
  AX_FCS_TYPE_REAL_SET([real*8],[MPI_REAL8],[${ax_check_kindof}],[${ax_kind_isoc}])
elif test "$1" = "real*16" ; then
  AX_CHECK_KINDOF([real*16])
  AX_FCS_TYPE_REAL_SET([real*16],[MPI_REAL16],[${ax_check_kindof}],[${ax_kind_isoc}])
else
  AC_MSG_ERROR([Datatype $1 is not an appropriate Fortran real type to use for FCS.])
fi
])


# Add configure options to select C and Fortran types.
AC_DEFUN_ONCE([AX_FCS_TYPES_ARGS],[
AC_ARG_ENABLE([fcs-int],
 [AS_HELP_STRING([--enable-fcs-int=TYPE],
  [set FCS C integer type @<:@int@:>@])])
AC_ARG_ENABLE([fcs-float],
 [AS_HELP_STRING([--enable-fcs-float=TYPE],
   [set FCS C floating type @<:@double@:>@])])
AC_ARG_ENABLE([fcs-integer],
 [AS_HELP_STRING([--enable-fcs-integer=TYPE],
  [set FCS Fortran integer type])])
AC_ARG_ENABLE([fcs-real],
 [AS_HELP_STRING([--enable-fcs-real=TYPE],
  [set FCS Fortran floating type])])
])


# Set up FCS integer and floating point types for C and Fortran.
AC_DEFUN([AX_FCS_TYPES],[

AC_REQUIRE([AX_FCS_TYPES_ARGS])

fcs_int=
fcs_integer=
case "x${enable_fcs_int}:x${enable_fcs_integer}" in
  x:x)
    fcs_int="int"
    ;;
  x*:x)
    fcs_int="${enable_fcs_int}"
    ;;
  x:x*)
    fcs_integer="${enable_fcs_integer}"
    ;;
  *)
    fcs_int="${enable_fcs_int}"
    fcs_integer="${enable_fcs_integer}"
    ;;
esac

fcs_float=
fcs_real=
case "x${enable_fcs_float}:x${enable_fcs_real}" in
  x:x)
    fcs_float="double"
    ;;
  x*:x)
    fcs_float="${enable_fcs_float}"
    ;;
  x:x*)
    fcs_real="${enable_fcs_real}"
    ;;
  *)
    fcs_float="${enable_fcs_float}"
    fcs_real="${enable_fcs_real}"
    ;;
esac

# Check if C and/or Fortran compilers are required to determine missing types.

# The following conditionals ensure that the AX_PROG_..._MPI macros determine only MPI C and/or Fortran compilers if they are required.
# If they are not required then the corresponding non-MPI C and Fortran compilers are set up by AX_PROG_..._MPI.

m4_define([have_f_not],[( test "x${enable_fcs_integer}" = x || test "x${enable_fcs_real}" = x )])
m4_define([have_c_not_or_have_f_not],[( test "x${enable_fcs_int}" = x || test "x${enable_fcs_float}" = x || test "x${enable_fcs_integer}" = x || test "x${enable_fcs_real}" = x )])
m4_define([have_c_not_and_have_f],[( ( test "x${enable_fcs_int}" = x && test "x${enable_fcs_integer}" != x ) || ( test "x${enable_fcs_float}" = x && test "x${enable_fcs_real}" != x ) )])

ax_fcs_types_sets_cc=
ax_fcs_types_sets_fc=

AC_PROVIDE_IFELSE([AC_PROG_FC],
 [AC_PROVIDE_IFELSE([AC_PROG_CC],
   [#echo "A: nothing"
   ],
   [AC_PROVIDE_IFELSE([AC_PROG_CXX],
     [#echo "B: require CC if have_c_not_or_have_f_not"
      ax_fcs_types_sets_cc=yes
      AC_REQUIRE([AX_PROG_CC_MPI],[AX_PROG_CC_MPI([have_c_not_or_have_f_not],,[have_c_not_or_have_f_not && AC_MSG_FAILURE([The FCS library requires an MPI C compiler.])])])
     ],
     [#echo "C: require CC if have_f_not"
      ax_fcs_types_sets_cc=yes
      AC_REQUIRE([AX_PROG_CC_MPI],[AX_PROG_CC_MPI([have_f_not],,[have_f_not && AC_MSG_FAILURE([The FCS library requires an MPI C compiler.])])])
     ]
    )
   ]
  )
 ],
 [AC_PROVIDE_IFELSE([AC_PROG_CC],
   [#echo "D: require FC if have_c_not_and_have_f"
    ax_fcs_types_sets_fc=yes
    AC_REQUIRE([AX_PROG_FC_MPI],[AX_PROG_FC_MPI([have_c_not_and_have_f],,[have_c_not_and_have_f && AC_MSG_FAILURE([The FCS library requires an MPI Fortran compiler.])])])
   ],
   [AC_PROVIDE_IFELSE([AC_PROG_CXX],
     [#echo "E: require CC and FC if have_c_not_and_have_f"
      ax_fcs_types_sets_cc=yes
      AC_REQUIRE([AX_PROG_CC_MPI],[AX_PROG_CC_MPI([have_c_not_and_have_f],,[have_c_not_and_have_f && AC_MSG_FAILURE([The FCS library requires an MPI C compiler.])])])
      ax_fcs_types_sets_fc=yes
      AC_REQUIRE([AX_PROG_FC_MPI],[AX_PROG_FC_MPI([have_c_not_and_have_f],,[have_c_not_and_have_f && AC_MSG_FAILURE([The FCS library requires an MPI Fortran compiler.])])])
     ],
     [#echo "F: nothing"
     ]
    )
   ]
  )
 ]
)

# Determine missing C and Fortran types.

# Don't do this if the C compiler was only set by us.
if test "x${ax_fcs_types_sets_cc}" != xyes ; then
  if test "x${fcs_int}" = x ; then :
    AX_FORTRAN_F2C_TYPE([fcs_int],[${fcs_integer}],[short,int,long,long long])
  fi
  if test "x${fcs_float}" = x ; then :
    AX_FORTRAN_F2C_TYPE([fcs_float],[${fcs_real}],[float,double,long double])
  fi
fi

# Don't do this if the Fortran compiler was only set by us.
if test "x${ax_fcs_types_sets_fc}" != xyes ; then
  if test "x${fcs_integer}" = x ; then :
    AX_FORTRAN_C2F_TYPE([fcs_integer],[${fcs_int}],[integer,integer*2,integer*4,integer*8])
  fi
  if test "x${fcs_real}" = x ; then :
    AX_FORTRAN_C2F_TYPE([fcs_real],[${fcs_float}],[real,real*4,real*8,real*16])
  fi
fi

# Set the FCS types.

# At this point AC_PROG_CC and AC_PROG_FC might only be provided, because we have used the corresponding _MPI macros.
# However, ax_fcs_types_sets_cc and ax_fcs_types_sets_fc tells us whether we have used the _MPI macros or not.
# Thus, we set FCS C types only if we have not used AC_PROG_CC_MPI and we set FCS Fortran types only if we have not used AC_PROG_FC_MPI.

ax_fcs_c_set=
if test "x${ax_fcs_types_sets_cc}" = x ; then
  AC_PROVIDE_IFELSE([AC_PROG_CC],
   [ax_fcs_c_set=yes
    AX_FCS_TYPE_INT([${fcs_int}])
    AX_FCS_TYPE_FLOAT([${fcs_float}])
    AC_MSG_NOTICE([FCS C types: ${fcs_int}, ${fcs_float}])
  ])
fi
if test "x${ax_fcs_c_set}" != xyes ; then :
  AC_PROVIDE_IFELSE([AC_PROG_CXX],
   [AX_FCS_TYPE_INT([${fcs_int}])
    AX_FCS_TYPE_FLOAT([${fcs_float}])
    AC_MSG_NOTICE([FCS C types: ${fcs_int}, ${fcs_float}])
  ])
fi
if test "x${ax_fcs_types_sets_fc}" = x ; then
  AC_PROVIDE_IFELSE([AC_PROG_FC],
   [AX_FCS_TYPE_INTEGER([${fcs_integer}],[${fcs_int}])
    AX_FCS_TYPE_REAL([${fcs_real}],[${fcs_float}])
    AC_MSG_NOTICE([FCS Fortran types: ${fcs_integer}, ${fcs_real}])
  ])
fi


# Modify configure arguments so that sub-configures do not need to determine missing types again.
if test "x${fcs_int}" != x ; then
  case ${ac_configure_args} in
    *--enable-fcs-int=*) ;;
    *) ac_configure_args="${ac_configure_args} '--enable-fcs-int=${fcs_int}'"
  esac
fi

if test "x${fcs_float}" != x ; then
  case ${ac_configure_args} in
    *--enable-fcs-float=*) ;;
    *) ac_configure_args="${ac_configure_args} '--enable-fcs-float=${fcs_float}'"
  esac
fi

if test "x${fcs_integer}" != x ; then
  case ${ac_configure_args} in
    *--enable-fcs-integer=*) ;;
    *) ac_configure_args="${ac_configure_args} '--enable-fcs-integer=${fcs_integer}'"
  esac
fi

if test "x${fcs_real}" != x ; then
  case ${ac_configure_args} in
    *--enable-fcs-real=*) ;;
    *) ac_configure_args="${ac_configure_args} '--enable-fcs-real=${fcs_real}'"
  esac
fi
])
