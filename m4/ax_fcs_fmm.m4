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

AC_DEFUN_ONCE([AX_FCS_FMM_ARGS],[

# Choose a communication library to use.
AC_ARG_ENABLE([fcs-fmm-comm],
 [AS_HELP_STRING([--enable-fcs-fmm-comm=COMM],
   [choose communication library COMM to use for FMM: a1 for A1 (only on BlueGene), armci for ARMCI, auto for automatic selection @<:@auto@:>@])],
 [],[enable_fcs_fmm_comm=auto])

# Limit the amount of unrolling done.
AC_ARG_ENABLE([fcs-fmm-unrolled],
  [AS_HELP_STRING([--enable-fcs-fmm-unrolled],
     [whether to build unrolled FMM operators @<:@no@:>@])],
  [], [enable_fcs_fmm_unrolled=no])

# Limit the maximum multipole moment.
m4_define([fmm_max_mpol_default], [30])
m4_define([fmm_max_mpol_min], [0])
m4_define([fmm_max_mpol_max], [50])
AC_ARG_ENABLE([fcs-fmm-max-mpol],
  [AS_HELP_STRING([--enable-fcs-fmm-max-mpol=N],
     [maximum FMM multipole moment ]fmm_max_mpol_min[<=N<=]fmm_max_mpol_max[
      @<:@]fmm_max_mpol_default[@:>@])],
  [if test "$enableval" = no; then
     enable_fcs_fmm_max_mpol=fmm_max_mpol_max
   fi],
  [enable_fcs_fmm_max_mpol=fmm_max_mpol_default])

])


AC_DEFUN([AX_FCS_FMM_TOP],[

AX_FCS_FMM_ARGS

#if $1 ; then
#  AC_MSG_NOTICE([do FMM solver stuff in top-level configure])
#fi

])


# Choose a communication library to use.
AC_DEFUN([AX_FCS_FMM_COMM],[
if test "x${enable_fcs_fmm_comm}" = xauto ; then
  case "${build_alias} ${host_alias} ${target_alias}" in
    *bgl*|*bgp*)
      enable_fcs_fmm_comm=a1
      ;;
    *bgq*)
      enable_fcs_fmm_comm=simple-armci
      ;;
    *)
      enable_fcs_fmm_comm=armci
      ;;
  esac
fi
case "${enable_fcs_fmm_comm}" in
  armci)
    AC_MSG_NOTICE([using ARMCI communication library])
    FMM_MP="FMM_MP_ARMCI"
    ;;
  a1)
    AC_MSG_NOTICE([using A1 communication library])
    FMM_MP="FMM_MP_A1"
    case $ac_configure_args in
      *--with-spi*) ;;
      *)
        ax_spi_includes="/bgsys/drivers/ppcfloor/arch/include"
        ax_spi_libs="/bgsys/drivers/ppcfloor/runtime/SPI"
        for ax_d in ${ax_spi_includes} ; do
          if test -d "${ax_d}" ; then
            AC_MSG_NOTICE([adding --with-spi-include=${ax_d} for A1 sub-configure])
            ac_configure_args="${ac_configure_args} --with-spi-include=${ax_d}"
            break
          fi
        done
        for ax_d in ${ax_spi_libs} ; do
          if test -d "${ax_d}" ; then
            AC_MSG_NOTICE([adding --with-spi-lib=${ax_d} for A1 sub-configure])
            ac_configure_args="${ac_configure_args} --with-spi-lib=${ax_d}"
            break
          fi
        done
        ;;
    esac
    ;;
  mpi)
    AC_MSG_FAILURE([MPI communication library cannot be used yet!])
    FMM_MP="FMM_MP_MPI"
    ;;
  simple-armci|sarmci)
    FMM_MP="FMM_MP_SIMPLE_ARMCI"
    case "${build_alias} ${host_alias} ${target_alias}" in
      *bgl*|*bgp*)
        enable_simple_armci_device=dcmfd
        ;;
      *bgq*)
        enable_simple_armci_device=pamid
        ;;
      *)
        enable_simple_armci_device=dcmfd
        ;;
    esac
    AC_MSG_NOTICE([using SIMPLE-ARMCI communication library with device '${enable_simple_armci_device}'])
    ;;
  *)
    AC_MSG_FAILURE([unknown communication library ${enable_fcs_fmm_comm} (use armci, a1 or auto)])
    ;;
esac
AC_DEFINE_UNQUOTED([FMM_MP],[${FMM_MP}],[Define to the communication library to use (FMM_MP_ARMCI or FMM_MP_A1).])
AM_CONDITIONAL(ENABLE_FMM_ARMCI,[test "x${FMM_MP}" = xFMM_MP_ARMCI -o "x${enable_dist}" = xyes])
AM_CONDITIONAL(ENABLE_FMM_A1,[test "x${FMM_MP}" = xFMM_MP_A1 -o "x${enable_dist}" = xyes])
AM_CONDITIONAL(ENABLE_FMM_MPI,[test "x${FMM_MP}" = xFMM_MP_MPI])
AM_CONDITIONAL(ENABLE_FMM_SIMPLE_ARMCI,[test "x${FMM_MP}" = xFMM_MP_SIMPLE_ARMCI -o "x${enable_dist}" = xyes])
if test "x${FMM_MP}" = xFMM_MP_ARMCI -o "x${enable_dist}" = xyes ; then
  AC_CONFIG_SUBDIRS([armci])
fi
if test "x${FMM_MP}" = xFMM_MP_A1 -o "x${enable_dist}" = xyes ; then
  AC_CONFIG_SUBDIRS([a1])
fi
AM_CONDITIONAL(ENABLE_SIMPLE_ARMCI_DCMFD,[test "x${enable_simple_armci_device}" = xdcmfd])
AM_CONDITIONAL(ENABLE_SIMPLE_ARMCI_PAMID,[test "x${enable_simple_armci_device}" = xpamid])
if test "x${FMM_MP}" = xFMM_MP_SIMPLE_ARMCI -o "x${enable_dist}" = xyes ; then
  AC_CONFIG_FILES([simple-armci/Makefile
    simple-armci/generic/Makefile
    simple-armci/dcmfd/Makefile
    simple-armci/pamid/Makefile])
fi
])


# Limit the amount of unrolling done.
AC_DEFUN([AX_FCS_FMM_UNROLLED],[
case $enable_fcs_fmm_unrolled in
yes)
  AC_MSG_NOTICE([enabling unrolled FMM operators])
  AC_DEFINE([FMM_UNROLLED], [],
    [enable unrolled FMM operators for compilation (significant increase in compile time)])
  ;;
*)
  AC_MSG_NOTICE([disabling unrolled FMM operators])
  ;;
esac
AM_CONDITIONAL(ENABLE_FMM_UNROLLED,[test "x$enable_fcs_fmm_unrolled" = xyes])
])


# Sanity checks: must be number, must be between allowed min and max.
AC_DEFUN([AX_FCS_FMM_MAX_MPOL],[
case $enable_fcs_fmm_max_mpol in
'' | *[[!0-9-]]*)
  AC_MSG_WARN([setting max multipole degree to default of fmm_max_mpol_default])
  enable_fcs_fmm_max_mpol=fmm_max_mpol_default ;;
*)
  if test $enable_fcs_fmm_max_mpol -lt 0; then
    AC_MSG_WARN([setting max multipole degree to minimum of fmm_max_mpol_min])
    enable_fcs_fmm_max_mpol=fmm_max_mpol_min
  elif test $enable_fcs_fmm_max_mpol -gt fmm_max_mpol_max; then
    AC_MSG_WARN([setting max multipole degree to maximum of fmm_max_mpol_max])
    enable_fcs_fmm_max_mpol=fmm_max_mpol_max
  else
    AC_MSG_NOTICE([setting max multipole degree to $enable_fcs_fmm_max_mpol])
  fi
esac
AC_DEFINE_UNQUOTED([FMM_MAXNMULTIPOLES], [$enable_fcs_fmm_max_mpol],
  [limit number of poles (implemented maximum=fmm_max_mpol_max])
])


# Test if Fortran 2003 procedure pointers are available
AC_DEFUN([AX_FCS_FMM_PROCPTR],[
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[[
      implicit none   
      procedure(test), pointer :: pptr => null()
      contains
      subroutine test()
      end subroutine test
]])],[],[disable_fmm_procedurepointers=yes])
AC_LANG_POP([Fortran])
case $disable_fmm_procedurepointers in
yes)
  AC_MSG_NOTICE([disabling usage of procedure pointers in FMM])
  AC_DEFINE([FMM_NOFUNCTIONPOINTER], [],
    [disable usage of procedure pointers in FMM])
  ;;
*)
  AC_MSG_NOTICE([enabling usage of procedure pointers in FMM])
  ;;
esac
])


AC_DEFUN([AX_FCS_FMM_SOLVER],[

AX_FCS_FMM_ARGS

AX_FCS_FMM_COMM

AX_FCS_FMM_UNROLLED

AX_FCS_FMM_MAX_MPOL

AX_FCS_FMM_PROCPTR

])
