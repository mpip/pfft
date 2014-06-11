
# AX_FORTRAN_DEPEND([DIR], [DEPFILE], [LOCAL_FILES], [FOREIGN_FILES], [OBJEXT], [MODEXI], [MOD_IGNORE], [MOD_TOUCH])
# ------------------------------------------------------------------------------------------------------------------

AC_DEFUN_ONCE([AX_FORTRAN_DEPEND_FUNC],[[
ax_fortran_depend_func() {

  local dir

  dir="${1}"

  local cur_builddir cur_srcdir sub_builddir sub_srcdir depfile

  cur_srcdir="`cd ${srcdir} ; pwd`"
  cur_builddir="`pwd`"
  sub_builddir="`cd ${dir} ; pwd`"
  sub_srcdir="`cd ${srcdir} ; cd ${dir} ; pwd`"

  depfile="${sub_builddir}/${2}"

  local loc_part_list for_part_list
  
  loc_part_list="${3}"
  for_part_list="${4}"

  local old_pwd UNIQSTR OBJEXT MODEXT SEDTMP FOREIGN

  old_pwd="`pwd`"

  UNIQSTR="QWERTZ"
  OBJEXT="${5:=o}"
  MODEXT="${6:=mod}"
  IGNORE="${7}"
  TOUCH="${8}"
  SEDTMP="${old_pwd}/qwertz_sed_tmp"
  FOREIGN=

  echo "ax_fortran_depend_func"
  echo "  dir: $dir"
  echo "  depfile: $depfile"
  echo "  cur_builddir: $cur_builddir"
  echo "  cur_srcdir: $cur_srcdir"
  echo "  sub_builddir: $sub_builddir"
  echo "  sub_srcdir: $sub_srcdir"
  echo "  loc_part_list: $loc_part_list"
  echo "  for_part_list: $for_part_list"
  echo "  UNIQSTR: $UNIQSTR"
  echo "  OBJEXT: $OBJEXT"
  echo "  MODEXT: $MODEXT"
  echo "  IGNORE: $IGNORE"
  echo "  TOUCH: $TOUCH"
  echo "  SEDTMP: $SEDTMP"
  echo "  FOREIGN: $FOREIGN"
  echo "  old_pwd: $old_pwd"

  ]AS_MKDIR_P($dir)[

  cd $sub_srcdir
  
  local f w n x d
  local uses useswo mods modswo uses_

  for x in ${IGNORE} ; do
    echo "s! ${x}\.${MODEXT}!!g"
  done > "${SEDTMP}"
  echo "/^.*:[ ]*$/d" >> "${SEDTMP}"

  echo '# This file is generated from configure.ac. Do not edit.
# Edit configure.ac instead. If the dependencies need to be updated,
# touch the configure.ac file and run make.
' > "${depfile}"

  for f in $loc_part_list ${UNIQSTR} ${for_part_list} ; do
#    echo "#f: ${f}"
    if test "${f}" = "${UNIQSTR}" ; then
#      echo "FOREIGN"
      FOREIGN=yes
      continue
    fi
    if test -f "${f}" ; then
#      echo "#f: ${f}"
#      cat ${f} 2>&1 | \
      cpp ${FCDEFS} ${FCCPPFLAGS} $CPPFLAGS -I${sub_builddir} -I${sub_srcdir} -I${cur_builddir} -I${cur_srcdir} ${f} 2>&1 | \
      sed -e "s/\!.*$//p" | \
      grep -e "use" -e "module" | \
      sed -n -e "s/^.*use[ \t][ \t]*\([a-zA-Z][a-zA-Z0-9_]*\).*$/${UNIQSTR}use \1/p" \
             -e "s/^.*module[ \t][ \t]*\([a-zA-Z][a-zA-Z0-9_]*\)$/${UNIQSTR}mod \1/p" \
             -e "\$i${UNIQSTR} ${UNIQSTR}" | \
      while read w n ; do
        if test "${w}" = "${UNIQSTR}use" ; then
          eval "useswo=\${uses##* ${n} *}"
          test "${uses}" = "${useswo}" && uses="${uses} ${n} "
        elif test "${w}" = "${UNIQSTR}mod" ; then
          eval "modswo=\${mods##* ${n} *}"
          test "${mods}" = "${modswo}" && mods="${mods} ${n} "
        elif test "${w} ${n}" = "${UNIQSTR} ${UNIQSTR}" ; then
          for x in ${uses} ; do
            eval "modswo=\${mods##* ${x} *}"
            test "${mods}" = "${modswo}" && uses_="${uses_} ${x} "
          done
          if test -n "${FOREIGN}" ; then
            d="${f%/*}"
            if test "${d}" = "${f}" ; then
              d="./"
            else
              d="${d}/"
            fi
          else
            d="./"
          fi
          if test -z "${FOREIGN}" ; then
            uses=${uses_}
            if test -n "${uses}" ; then
              printf "${f%.*}.${OBJEXT}:"
              for x in ${uses} ; do
                printf " ${x}.${MODEXT}"
              done
              echo ""
            fi
            for x in ${mods} ; do
              echo "${x}.${MODEXT}: ${f%.*}.${OBJEXT}"
              case " ${TOUCH} " in
                *\ ${x}\ *)
                echo "	touch $""@"
                ;;
              esac
            done
          fi
          for x in ${mods} ; do
            test -n "${d#./}" && echo "s!${x}.${MODEXT}!${d#./}${x}.${MODEXT}!g" >> "${SEDTMP}"
          done
        fi
      done
    fi
  done >> ${depfile}

  sed -i -f "${SEDTMP}" "${depfile}"

  rm -f "${SEDTMP}"

  cd $old_pwd
}
]])

AC_DEFUN([AX_FORTRAN_DEPEND],[
AC_PROG_CPP
AX_FORTRAN_DEPEND_FUNC
ax_fortran_depend_func "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
])


# AX_FORTRAN_MOD_DEPEND([FILE], [DIR], [LOCAL_FILES], [FOREIGN_FILES], [OBJEXT], [MODEXI], [MOD_IGNORE], [MOD_TOUCH])
# -------------------------------------------------------------------------------------------------------------------

AC_DEFUN([AX_FORTRAN_MOD_DEPEND],[
  m4_define([loc_pat_list],m4_bpatsubst($3,[[*]+],[\\*]))
  m4_define([for_pat_list],m4_bpatsubst($4,[[*]+],[\\*]))
  m4_pushdef([fortran_mod_depend_build_error],[m4_esyscmd([
ax_UNIQSTR="QWERTZ"
ax_OBJEXT=]m4_ifblank($5,[o],$5)[
ax_MODEXT=]m4_ifblank($6,[mod],$6)[
ax_SEDTMP="qwertz_sed_tmp"
ax_FOREIGN=
cd $2 && echo "" > "${ax_SEDTMP}" &&
echo '# This file is generated from configure.ac. Do not edit.
# Edit configure.ac instead. If the dependencies need to be updated,
# touch the configure.ac file and run make.
' > $1 && 
for ax_p in ]loc_pat_list[ ${ax_UNIQSTR} ]for_pat_list[ ; do
  if test "${ax_p}" = "${ax_UNIQSTR}" ; then
#    echo "FOREIGN"
    ax_FOREIGN=yes
    continue
  fi
#  echo "#p: ${ax_p}"
  eval ax_f_list="${ax_p}"
  for ax_f in ${ax_f_list} ; do
    if test -f "${ax_f}" ; then
#      echo "#f: ${ax_f}"
      cat ${ax_f} | \
      sed -e "s/\!.*$//p" | \
      grep -e "use" -e "module" | \
      sed -n -e "s/^.*use[ \t][ \t]*\([a-zA-Z][a-zA-Z0-9_]*\).*$/${ax_UNIQSTR}use \1/p" \
             -e "s/^.*module[ \t][ \t]*\([a-zA-Z][a-zA-Z0-9_]*\)$/${ax_UNIQSTR}mod \1/p" \
             -e "\$i${ax_UNIQSTR} ${ax_UNIQSTR}" | \
      while read ax_w ax_n ; do
        if test "${ax_w}" = "${ax_UNIQSTR}use" ; then
          eval "ax_useswo=\${ax_uses##* ${ax_n} *}"
          test "${ax_uses}" = "${ax_useswo}" && ax_uses="${ax_uses} ${ax_n} "
        elif test "${ax_w}" = "${ax_UNIQSTR}mod" ; then
          eval "ax_modswo=\${ax_mods##* ${ax_n} *}"
          test "${ax_mods}" = "${ax_modswo}" && ax_mods="${ax_mods} ${ax_n} "
        elif test "${ax_w} ${ax_n}" = "${ax_UNIQSTR} ${ax_UNIQSTR}" ; then
          for ax_x in ${ax_uses} ; do
            eval "ax_modswo=\${ax_mods##* ${ax_x} *}"
            test "${ax_mods}" = "${ax_modswo}" && ax_uses_="${ax_uses_} ${ax_x} "
          done
          if test -n "${ax_FOREIGN}" ; then
            ax_d="${ax_f%/*}"
            if test "${ax_d}" = "${ax_f}" ; then
              ax_d="./"
            else
              ax_d="${ax_d}/"
            fi
          else
            ax_d="./"
          fi
          if test -z "${ax_FOREIGN}" ; then
            ax_uses=${ax_uses_}
            if test -n "${ax_uses}" ; then
              printf "${ax_f%.*}.${ax_OBJEXT}:"
              for ax_x in ${ax_uses} ; do
                printf " ${ax_x}.${ax_MODEXT}"
              done
              echo ""
            fi
            for ax_x in ${ax_mods} ; do
              echo "${ax_x}.${ax_MODEXT}: ${ax_f%.*}.${ax_OBJEXT}"
              case " $8 " in
                *\ ${ax_x}\ *)
                echo "	touch $""@"
                ;;
              esac
            done
          fi
          for ax_x in ${ax_mods} ; do
            test -n "${ax_d#./}" && echo "s!${ax_x}.${ax_MODEXT}!${ax_d#./}${ax_x}.${ax_MODEXT}!g" >> "${ax_SEDTMP}"
          done
        fi
      done
    fi
  done
done | sed]m4_foreach_w([ignore],[$7],[ -e "s! ]ignore[\.${ax_MODEXT}!!g"])[ -e "/^.*:[ ]*$/d" >> $1 && sed -i -f "${ax_SEDTMP}" $1 && rm -f "${ax_SEDTMP}"])])dnl
  m4_ifval(fortran_mod_depend_build_error,
    [m4_fatal([generating Fortran module dependencies failed: ]fortran_mod_depend_build_error[.])])dnl
  m4_popdef([fortran_mod_depend_build_error])dnl
])
