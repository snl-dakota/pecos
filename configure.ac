#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.60])
AC_INIT([PECOS],[1.0],[pecos-developers@development.sandia.gov])

AC_CONFIG_SRCDIR([src/pecos_data_types.hpp])
AC_CONFIG_HEADER([pecos_config.h])

AM_INIT_AUTOMAKE([1.9.6 foreign tar-ustar])

AC_LANG([C++])

# Future feature (for now, we rely on ACLOCAL_AMFLAGS in Makefile.am):
#AC_CONFIG_MACRO_DIR([../../m4])

# Since we're using libtool, don't build shared libraries!
dnl WJB: ASAP -- need to do some more build testing -- AC_DISABLE_SHARED

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_INSTALL

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL

# Checks for libraries.
AC_LANG_PUSH([C])
#### WJB: really needed?? AC_CHECK_LIB([m])
AC_LANG_POP([C])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([cstdlib])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
## WJB: later -- AC_CHECK_FUNCS([rand random])


# PECOS-specific package checks.
RAC_PATH_GSL
RAC_PATH_TEUCHOS


# Write the Makefiles and exit.
AC_CONFIG_FILES([Makefile packages/Makefile src/Makefile])

AC_OUTPUT
