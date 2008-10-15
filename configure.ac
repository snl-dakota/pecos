#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.60])
AC_INIT([PECOS],[1.0],[pecos-developers@development.sandia.gov])
AM_INIT_AUTOMAKE([1.9.6 foreign tar-ustar])
AC_LANG([C++])

AC_CONFIG_SRCDIR([src/pecos_data_types.hpp])
AC_CONFIG_HEADER([pecos_config.h])

# Future feature (for now, we rely on ACLOCAL_AMFLAGS in Makefile.am):
#AC_CONFIG_MACRO_DIR([m4])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
dnl AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_ARG_ENABLE([f90],AS_HELP_STRING([--disable-f90],[turn F90 support off]),
	      [enable_f90=$enableval],[enable_f90=yes])
AM_CONDITIONAL([ENABLE_F90],[test "x$enable_f90" = xyes])
# Checks for typedefs, structures, and compiler characteristics.
AC_F77_WRAPPERS
if test "x$enable_f90" = xyes; then
  AC_PROG_FC([],[90])
  AC_FC_WRAPPERS
  AC_LANG_PUSH([Fortran])
  AC_FC_SRCEXT([f90])
  AC_FC_FIXEDFORM
  AC_LANG_POP([Fortran])
fi
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S

# Pecos is a C++ code.
AC_LANG([C++])

# Checks for header files.
AC_LANG_PUSH([C])
AC_HEADER_STDC
AC_LANG_POP([C])
AC_CHECK_HEADERS([cstdlib])
AC_HEADER_STDBOOL


ACX_BLAS
ACX_LAPACK

# Checks for library functions.
## WJB: later -- AC_CHECK_FUNCS([rand random])


# PECOS-specific package checks.
PECOS_PACKAGES


# Write the Makefiles and exit.
AC_CONFIG_FILES([Makefile packages/Makefile src/Makefile test/Makefile])

AC_OUTPUT
