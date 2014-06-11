
AC_DEFUN([AX_SORT_LIB],[

# Get the MPI C compiler.
AC_REQUIRE([AX_PROG_CC_MPI],[AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The Sorting Library requires an MPI C compiler.]))])

# Try to put the C compiler in C99 mode (this macro is relatively new in Autoconf).
#m4_ifdef([AC_PROG_CC_C99],
#  [AC_REQUIRE([AC_PROG_CC_C99])])

# Redefine "inline" if the C compiler has problems with it.
AC_REQUIRE([AC_C_INLINE])

# We need the math library.
AC_SEARCH_LIBS([log], [m])

# Blue Gene specific checks
AC_CHECK_HEADERS([spi/kernel_interface.h common/bgp_personality.h common/bgp_personality_inlines.h],
  [], [],
  [[#ifdef HAVE_SPI_KERNEL_INTERFACES_H
    # include <spi/kernel_interface.h>
    #endif
    #ifdef HAVE_COMMON_BGP_PERSONALITY_H
    # include <common/bgp_personality.h>
    #endif]])

AC_CHECK_TYPES([_BGP_Personality_t],
  [], [],
  [AC_INCLUDES_DEFAULT
   [#ifdef HAVE_SPI_KERNEL_INTERFACES_H
    # include <spi/kernel_interface.h>
    #endif
    #ifdef HAVE_COMMON_BGP_PERSONALITY_H
    # include <common/bgp_personality.h>
    #endif
    #ifdef HAVE_COMMON_BGP_PERSONALITY_INLINES_H
    # include <common/bgp_personality_inlines.h>
    #endif]])

])
