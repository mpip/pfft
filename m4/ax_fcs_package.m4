# Copyright (C) 2011 The ScaFaCoS project
#  
# This file is part of ScaFaCoS.
#  
# ScaFaCoS is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
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

AC_DEFUN([AX_FCS_PACKAGE_RESET],[
  ax_fcs_package_file="fcs-package.info"
  test -n "$1" && ax_fcs_package_file="$1/${ax_fcs_package_file}"
  rm -f ${ax_fcs_package_file}
])


AC_DEFUN([AX_FCS_PACKAGE_ADD],[
  ax_fcs_package_file="fcs-package.info"
  test -n "$3" && ax_fcs_package_file="$3/${ax_fcs_package_file}"
  echo "$1=$2" >> ${ax_fcs_package_file}
])


AC_DEFUN([AX_FCS_PACKAGE_PROCESS],[
  ax_fcs_topdir="$1"
  ax_fcs_pkgdirs="$2"
  ax_fcs_package_file="fcs-package.info"
  ax_nl="
"
  test "x${ax_fcs_topdir}" != x && ax_fcs_topdir="${ax_fcs_topdir%/}/"
  while test -n "${ax_fcs_pkgdirs}" ; do
#    echo "packages: ${ax_fcs_pkgdirs}"
    ax_p="${ax_fcs_pkgdirs%% *}"
    if test "${ax_p}" = "${ax_fcs_pkgdirs}" ; then
      ax_fcs_pkgdirs=
    else
      ax_fcs_pkgdirs="${ax_fcs_pkgdirs#* }"
    fi
    ax_p="${ax_p%/}/"
    ax_p="${ax_p%./}"
#    echo "ax_p: $ax_p"
    ax_f="${ax_fcs_topdir}${ax_p}${ax_fcs_package_file}"
#    echo "ax_f: ${ax_f}"
    if test -f "${ax_f}" ; then
      ax_cfg=`cat ${ax_f}`
#      echo "$ax_cfg"
      _ifs="${IFS}"
      IFS="${ax_nl}"
      for ax_l in $ax_cfg ; do
        IFS="${_ifs}"
#        echo "ax_l: ${ax_l}"
        ax_var="${ax_l%%=*}"
        ax_val="${ax_l#*=}"
#        echo "-> $ax_var / $ax_val"
        test "x${ax_var}" = x -o "x${ax_val}" = x && continue
        if test "x${ax_var}" = xSUB_PACKAGES ; then
          for ax_x in ${ax_val} ; do
            ax_fcs_pkgdirs="${ax_p}${ax_x} ${ax_fcs_pkgdirs}"
          done
          continue
        fi
        if test "x${ax_var%LIBS_A}" != "x${ax_var}" ; then
          ax_val_new=
          for ax_x in ${ax_val} ; do
            if test "x${ax_x#/}" = "x${ax_x}" ; then
              ax_x="${ax_p}${ax_x}"
            fi
            ax_val_new="${ax_val_new} ${ax_x}"
          done
          ax_val="${ax_val_new# }"
        fi
        eval ax_old=\"\$ax_fcs_package_${ax_var}\"
#        echo "-> old: ${ax_old}"
        if test "x${ax_old}" = x ; then
          eval ax_fcs_package_${ax_var}=\"${ax_val}\"
        else
          eval ax_fcs_package_${ax_var}=\"\$ax_fcs_package_${ax_var} ${ax_val}\"
        fi
#        eval ax_new=\"\$ax_fcs_package_${ax_var}\"
#        echo "-> new: ${ax_new}"
        IFS="${ax_nl}"
      done
      IFS="${_ifs}"
    fi
  done
])
