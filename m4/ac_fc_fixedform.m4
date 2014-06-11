# This macro has been submitted for inclusion into the Autoconf 2.66 tree.
# Use the upstream version when it is defined.
m4_version_prereq([2.66], [],
[

# AC_FC_FIXEDFORM([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept
# fixed-format source code, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using new extension) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
#
# The known flags are:
#       -ffixed-form: GNU g77, gfortran, g95
#             -fixed: Intel compiler (ifort), Sun compiler (f95)
#            -qfixed: IBM compiler (xlf*)
#            -Mfixed: Portland Group compiler
#         -fixedform: SGI compiler
#           -f fixed: Absoft Fortran
#      +source=fixed: HP Fortran
#    (-)-fix, -Fixed: Lahey/Fujitsu Fortran
#             -fixed: NAGWare
# Since compilers may accept fixed form based on file name extension,
# but users may want to use it with others as well, call AC_FC_SRCEXT
# with the respective source extension before calling this macro.
AC_DEFUN_ONCE([AC_FC_FIXEDFORM],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK([for Fortran flag needed to accept fixed-form source],
	       [ac_cv_fc_fixedform],
[ac_cv_fc_fixedform=unknown
ac_fc_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffixed-form -fixed -qfixed -Mfixed -fixedform "-f fixed" \
	       +source=fixed -fix --fix -Fixed
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_fixedform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[
C     This comment should confuse free-form compilers.
      program main
      end]],
		    [ac_cv_fc_fixedform=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_fixedform_FCFLAGS_save
])
if test "x$ac_cv_fc_fixedform" = xunknown; then
  m4_default([$2],
	     [AC_MSG_ERROR([Fortran does not accept fixed-form source], 77)])
else
  if test "x$ac_cv_fc_fixedform" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_fixedform"
  fi
  $1
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_FIXEDFORM

])dnl m4_version_prereq([2.66])
