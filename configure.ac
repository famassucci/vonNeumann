#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(vonNeumann, 0.1.0, francesco.massucci@urv.cat)
AC_CONFIG_SRCDIR([src/optimal_flux.c])
AC_CONFIG_HEADERS([config.h])

PACKAGE=vonNeumann
VERSION=0.1.0

# Checks for programs.
AC_PROG_CC
AM_INIT_AUTOMAKE
LT_INIT
AM_PROG_CC_C_O

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

# Checks for libraries.
AC_CHECK_LIB(c, main)
AC_CHECK_LIB(m, [sqrt])

# Global CFLAGS
WARN_FLAGS="-Wall -Wextra -Wshadow -Wno-variadic-macros --pedantic"
CFLAGS="$WARN_FLAGS -O2"
CFLAGS_save="$CFLAGS"

# enable "optimize option", which will compile using only the -O3 option
AC_ARG_ENABLE([optimize],
    [AC_HELP_STRING([--enable-optimize],
        [disable debugging (removes -g flags), disable compiler warnings\
	     and set optimization level O3 (disabled)])],
    [CFLAGS="-O3"], [CFLAGS="$CFLAGS_save"])

# enable "more-debug", which will compile using extra debugging options
CFLAGS_save="$CFLAGS"
AC_ARG_ENABLE([more-debug],
    [AC_HELP_STRING([--enable-more-debug],
        [enable more debugging information (disabled)])],
     [CFLAGS="$WARN_FLAGS -g3 -ggdb -O0"], [CFLAGS="$CFLAGS_save"])

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h math.h])
VN_CFLAGS="-I$includedir"
VN_LIBS="-L$libdir -lvonNeumann"

AC_SUBST(VN_CFLAGS)
AC_SUBST(VN_LIBS)


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

# output
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 vonNeumann/Makefile
		 vonNeumann-config vonNeumann.pc])
AC_OUTPUT
