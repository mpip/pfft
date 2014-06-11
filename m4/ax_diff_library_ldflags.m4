
AC_DEFUN([AX_DIFF_LIBRARY_LDFLAGS_CXX],[
AX_CC_LIBRARY_LDFLAGS
AX_CXX_LIBRARY_LDFLAGS
AC_CACHE_CHECK([difference between LDFLAGS of C++ and C.],ax_cv_diff_library_ldflags_cxx,[
ax_clibs=$CLIBS
ax_cxxlibs=" $CXXLIBS "
#echo "ax_clibs: $ax_clibs"
#echo "ax_cxxlibs: $ax_cxxlibs"
for ax_l in ${ax_clibs} ; do
  ax_cxxlibs=`echo "${ax_cxxlibs}" | sed -e "s! ${ax_l} ! !g"`
done
ax_cv_diff_library_ldflags_cxx=
for ax_x in ${ax_cxxlibs} ; do ax_cv_diff_library_ldflags_cxx="${ax_cv_diff_library_ldflags_cxx} ${ax_x}" ; done
ax_cv_diff_library_ldflags_cxx="${ax_cv_diff_library_ldflags_cxx# }"
])
CXXLIBS_DIFF=${ax_cv_diff_library_ldflags_cxx}
AC_SUBST([CXXLIBS_DIFF])
])


AC_DEFUN([AX_DIFF_LIBRARY_LDFLAGS_FC],[
AX_CC_LIBRARY_LDFLAGS
AC_FC_LIBRARY_LDFLAGS
m4_ifdef([ACX_FC_LIBRARY_LDFLAGS_FIX],[ACX_FC_LIBRARY_LDFLAGS_FIX])
AC_CACHE_CHECK([difference between LDFLAGS of Fortran and C.],ax_cv_diff_library_ldflags_fc,[
ax_clibs=$CLIBS
ax_fclibs=" $FCLIBS "
#echo "ax_clibs: $ax_clibs"
#echo "ax_fclibs: $ax_fclibs"
for ax_l in ${ax_clibs} ; do
  ax_fclibs=`echo "${ax_fclibs}" | sed -e "s! ${ax_l} ! !g"`
done
ax_cv_diff_library_ldflags_fc=
for ax_x in ${ax_fclibs} ; do ax_cv_diff_library_ldflags_fc="${ax_cv_diff_library_ldflags_fc} ${ax_x}" ; done
ax_cv_diff_library_ldflags_fc="${ax_cv_diff_library_ldflags_fc# }"
])
FCLIBS_DIFF=${ax_cv_diff_library_ldflags_fc}
AC_SUBST([FCLIBS_DIFF])
])
