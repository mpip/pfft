# ACX_MPIEXEC([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
# A few heuristics to find out how to start an MPI job.
# This is currently very incomplete, see test/defs.in for
# helper functions which do the actual work for the known
# set of systems.
#
# We should probably have facilities to specify not only
# the total number of processes but also the number of
# processes per node etc.
AC_DEFUN([ACX_MPIEXEC], [
AC_REQUIRE([ACX_MPI])dnl

AC_CHECK_PROGS([MPIEXEC], [llrun msub qsub mpiexec mpirun mpprun], [false])
AC_ARG_VAR([MPIEXEC], [program to start MPI jobs on the host system])

AC_CACHE_CHECK([how to pass the number of processes to $MPIEXEC],
	       [acx_cv_mpiexec_minus_n], [
if test "${MPIEXEC_MINUS_N+set}" = set; then
  acx_cv_mpiexec_minus_n=$MPIEXEC_MINUS_N
else
  case $MPIEXEC in #(
    mpiexec|*/mpiexec)  acx_cv_mpiexec_minus_n="-n " ;; #(
    mpprun|*/mpprun)    acx_cv_mpiexec_minus_n="-n " ;; #(
    mpirun|*/mpirun)    acx_cv_mpiexec_minus_n="-np " ;; #(
    llrun|*/llrun)      acx_cv_mpiexec_minus_n="-np ";; #(
    msub|*/msub)        acx_cv_mpiexec_minus_n="" ;; # write a script instead (
    qsub|*/qsub)        acx_cv_mpiexec_minus_n="" ;; # write a script instead (
    *)                  acx_cv_mpiexec_minus_n=unknown ;;
  esac
fi
])
if test "x$acx_cv_mpiexec_minus_n" = xunknown; then
  MPIEXEC_MINUS_N=
else
  MPIEXEC_MINUS_N=$acx_cv_mpiexec_minus_n
fi
AC_SUBST([MPIEXEC_MINUS_N])dnl
])dnl ACX_MPIEXEC
