# based on: http://www.ccp4.ac.uk/dist/x-windows/Mosflm/macros/cxx.m4

AC_DEFUN([_AX_PROG_CXX_V_OUTPUT],
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_LANG_PUSH(C++)dnl
AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
ax_save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS m4_default([$1], [$ax_cv_prog_cxx_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ax_cxx_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ax_cxx_v_output" >&AS_MESSAGE_LOG_FD
CXXFLAGS=$ax_save_CXXFLAGS

rm -f conftest*
AC_LANG_POP(C++)dnl    
# If we are using xlC then replace all the commas with spaces.    
if echo $ax_cxx_v_output | grep xlC >/dev/null 2>&1; then
  ax_cxx_v_output=`echo $ax_cxx_v_output | sed 's/,/ /g'`
fi
# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ax_cxx_v_output="`echo $ax_cxx_v_output |
        grep 'LPATH is:' |
        sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ax_cxx_v_output"

])# _AX_PROG_CXX_V_OUTPUT


AC_DEFUN([_AX_PROG_CXX_V],
[AC_CACHE_CHECK([how to get verbose linking output from $CXX],
                [ax_cv_prog_cxx_v],
[AC_LANG_ASSERT(C++)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ax_cv_prog_cxx_v=
# Try some options frequently used verbose output
for ax_verb in -v -verbose --verbose -V -\#\#\#; do
  _AX_PROG_CXX_V_OUTPUT($ax_verb)
  # look for -l* and *.a constructs in the output
  for ax_arg in $ax_cxx_v_output; do
     case $ax_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ax_cv_prog_cxx_v=$ax_verb
          break 2 ;;
     esac
  done
done
if test -z "$ax_cv_prog_cxx_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $CXX])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AX_PROG_CXX_V


AC_DEFUN([AX_CXX_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(C++)dnl
_AX_PROG_CXX_V  
AC_CACHE_CHECK([for C++ libraries], ax_cv_cxxlibs,
[if test "x$CXXLIBS" != "x"; then
  ax_cv_cxxlibs="$CXXLIBS" # Let the user override the test.
else

_AX_PROG_CXX_V_OUTPUT     

# If we find a line that looks like a linker call "ld ...", restrict analysis of linker options to that line.
ax_x=`echo "$ax_cxx_v_output" | grep "^ld "`
if test "x${ax_x}" != x ; then
  ax_cxx_v_output="${ax_x}"
fi

ax_cv_cxxlibs=

# Save positional arguments (if any)  
ax_save_positional="$[@]"

set X $ax_cxx_v_output
while test $[@%:@] != 1; do      
  shift
  ax_arg=$[1]
  case $ax_arg in   
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ax_arg, $ax_cv_cxxlibs, ,
              ax_cv_cxxlibs="$ax_cv_cxxlibs $ax_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ax_arg, $ax_cv_cxxlibs, ,
             [_AC_LINKER_OPTION([$ax_arg], ax_cv_cxxlibs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lcrt1.o | -lcrt2.o |-lcrtbegin.o )
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ax_cv_cxxlibs="$ax_cv_cxxlibs $ax_arg"
          ;;
        -[[LRuY]])      
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
          case $[2] in
             "" | -*);;       
             *)
                ax_arg="$ax_arg$[2]"
                shift; shift
                set X $ax_arg "$[@]"
                ;;
          esac
          ;;
        -YP,*)
          for ax_j in `echo $ax_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ax_j, $ax_cv_cxxlibs, ,
                               [ax_arg="$ax_arg $ax_j"
                               ax_cv_cxxlibs="$ax_cv_cxxlibs $ax_j"]) 
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ax_arg, $ax_cv_cxxlibs, ,
                             ax_cv_cxxlibs="$ax_cv_cxxlibs $ax_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments    
set X $ax_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ax_ld_run_path=`echo $ax_cxx_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ax_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ax_ld_run_path], ax_cv_cxxlibs)
      ;;
esac
fi # test "x$CXXLIBS" = "x"
])
CXXLIBS="$ax_cv_cxxlibs"             
AC_SUBST(CXXLIBS)
AC_LANG_POP(C++)dnl
])# AX_CXX_LIBRARY_LDFLAGS
