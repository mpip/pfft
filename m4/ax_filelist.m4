
# AX_FILELIST(FILE, [DIR], PATTERNS, [IGNORE-PATTERNS], [VAR], [IF-COND])
# -----------------------------------------------------------------------
# Create file FILE in directory DIR containing a list with all files that
# correspond to the file patterns in PATTERN (with respect to directory DIR).
# Files that correspond to the file patterns in IGNORE-PATTERNS are ignored.
# The file list variable is named VAR (or 'filelist' if VAR is empty).
# If IF-COND is not empty, then each file is added separately to the list using
# the conditional statement 'if IF-COND ... endif'. The current file name and 
# the number of the file is available in variables ax_filelist_name and
# ax_filelist_num, respectively.

AC_DEFUN([AX_FILELIST],[
  m4_define([pat_list],m4_bpatsubst($3,[[*]+],[\\*]))
  m4_define([ign_list],m4_bpatsubst($4,[[*]+],[\\*]))
  m4_define([var_name],m4_ifblank($5,filelist,$5))
  m4_pushdef([filelist_build_error],[m4_esyscmd([
cd ]m4_ifblank($2,./,$2)[ &&
echo '# This file is generated from configure.ac. Do not edit.
# Edit configure.ac instead. If the list of files needs to be updated,
# touch the configure.ac file and run make.
' > $1 &&
]m4_ifblank($6,[printf],[echo]) var_name[' =' >> $1 &&
for ax_p in ]m4_ifblank(ign_list,"",ign_list)[ ; do
#  echo "p: ${ax_p}"
  eval ax_f_list="${ax_p}"
  for ax_f in ${ax_f_list} ; do
#    echo "f: ${ax_f}"
    ax_f_list_ign="${ax_f_list_ign} ${ax_f}"
  done
done >> $1 &&
ax_filelist_num=0 &&
for ax_p in ]m4_ifblank(pat_list,"",pat_list)[ ; do
#  echo "p: ${ax_p}"
  eval ax_f_list="${ax_p}"
  for ax_filelist_name in ${ax_f_list} ; do
    case " ${ax_f_list_ign} " in
      *\ ${ax_filelist_name}\ *)
#        echo "ignoring ${ax_filelist_name}"
        ;;
      *)
        if test -f "${ax_filelist_name}" ; then
          ax_filelist_num=$((ax_filelist_num+1))
#          echo "f: ${ax_filelist_name}"
          ]m4_ifblank($6,[
          echo " \\"
          printf "${ax_filelist_name}"
            ],[
          echo "if $6"
          echo " var_name += ${ax_filelist_name}"
          echo "endif"
            ])[
        fi
        ;;
    esac
  done
done >> $1])])dnl
  m4_ifval(filelist_build_error,
    [m4_fatal([generating filelists failed: ]filelist_build_error[.])])dnl
  m4_popdef([filelist_build_error])dnl
])
