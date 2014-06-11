
AC_DEFUN([AX_CHECK_IBM_INTRINSICS],[
AC_CACHE_CHECK(
 [for IBM intrinsics in Fortran],
 [ax_cv_have_ibm_f_intrinsics],
 [save_FCFLAGS=$FCFLAGS
  FCFLAGS="$FCFLAGS -qarch=450d"
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([],[[
      real(kind=8) :: scr(2)
      complex(kind=8) :: reg
      call alignx(16, scr)
      reg = loadfp(scr(1))
]])],[ax_cv_have_ibm_f_intrinsics=yes],[ax_cv_have_ibm_f_intrinsics=no])
  AC_LANG_POP([Fortran])
  FCFLAGS=$save_FCFLAGS
])
if test "x${ax_cv_have_ibm_f_intrinsics}" = xyes ; then
  AC_DEFINE([HAVE_IBM_F_INTRINSICS],[1],[Define if IBM intrinsics can be used in Fortran.])
fi
ax_intrinsics_ibm_f="${ax_cv_have_ibm_f_intrinsics}"
])

AC_DEFUN([AX_CHECK_SSE_INTRINSICS],[
AC_CACHE_CHECK(
 [for SSE intrinsics in C],
 [ax_cv_have_sse_c_intrinsics],
 [AC_LINK_IFELSE([
  AC_LANG_PROGRAM([[
#include <xmmintrin.h>
__m128 testfunc(float *a, float *b) {
return _mm_add_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
}
]])],[ax_cv_have_sse_c_intrinsics=yes],[ax_cv_have_sse_c_intrinsics=no])
])
if test "x${ax_cv_have_sse_c_intrinsics}" = xyes ; then
  AC_DEFINE([HAVE_SSE_C_INTRINSICS],[1],[Define if SSE intrinsics can be used in C.])
fi
ax_intrinsics_sse_c="${ax_cv_have_sse_c_intrinsics}"
])

AC_DEFUN([AX_CHECK_SSE2_INTRINSICS],[
AC_CACHE_CHECK(
 [for SSE2 intrinsics in C],
 [ax_cv_have_sse2_c_intrinsics],
 [AC_LINK_IFELSE([
  AC_LANG_PROGRAM([[
#include <emmintrin.h>
__m128d testfunc(double *a, double *b) {
return _mm_add_pd(_mm_loadu_pd(a), _mm_loadu_pd(b));
}
]])],[ax_cv_have_sse2_c_intrinsics=yes],[ax_cv_have_sse2_c_intrinsics=no])
])
if test "x${ax_cv_have_sse2_c_intrinsics}" = xyes ; then
  AC_DEFINE([HAVE_SSE2_C_INTRINSICS],[1],[Define if SSE2 intrinsics can be used in C.])
fi
ax_intrinsics_sse2_c="${ax_cv_have_sse2_c_intrinsics}"
])
