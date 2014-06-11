
AC_DEFUN([AX_DISABLE],[
  m4_foreach([id], [$1], [
    m4_define([id_esc],AS_TR_SH(id))
    m4_define([id_var],ax_disable_[]id_esc)
    id_var=no
    AC_ARG_ENABLE(id,AS_HELP_STRING([--disable-]id,[disable ]id),
      [test "x$enableval" = xno && id_var=yes], [])
#    echo "id_var: $id_var"
  ])
])
