# ACX_FC_LIBRARY_LDFLAGS_FIX
# --------------------------
# Small fix to AC_FC_LIBRARY_LDFLAGS needed for xlf2003_r on BlueGene.
AC_DEFUN([ACX_FC_LIBRARY_LDFLAGS_FIX],
[new_FCLIBS=
for flag in $FCLIBS
do
  case $flag in #(
  -link) ;; # ignore this.  (
  *) new_FCLIBS="$new_FCLIBS $flag" ;;
  esac
done
FCLIBS=$new_FCLIBS
])
