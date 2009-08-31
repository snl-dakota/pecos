dnl PECOS Options

AC_DEFUN([PECOS_OPTIONS],[
  dnl Debug option check.
  AC_ARG_ENABLE([debug],
                AS_HELP_STRING([--enable-debug],[turn debug support on]),
                [enable_debug=$enableval],[enable_debug=no])
  dnl AC_MSG_CHECKING([debug])
  if test "x$enable_debug" = xyes; then
    dnl AC_MSG_RESULT([yes])
    dnl AC_DEFINE([HAVE_PECOS_DEBUG],[1],[Macro to enable debug in Pecos.])
    DEBUGFLAGS="-g"
  else
    dnl AC_MSG_RESULT([no])
    AC_DEFINE([BOOST_DISABLE_ASSERTS],[1],[Macro to disable Boost asserts.])
    DEBUGFLAGS=""
  fi
  AC_SUBST(DEBUGFLAGS)
  AM_CONDITIONAL([WITH_DEBUG_ENABLED],[test "x$enable_debug" = xyes])


  dnl Specify optimization level.
  AC_ARG_WITH([opt],
              AS_HELP_STRING([--with-opt=NUM],[set opt level (2 by default)]),
              [OPT_LEVEL=$withval],[OPT_LEVEL=2])
  AC_MSG_CHECKING([optimization level])
  if test "x$with_opt" = xno; then
    AC_MSG_RESULT([0])
    OPTFLAGS=""
  elif test "x$with_opt" = xyes; then
    AC_MSG_RESULT([2])
    OPTFLAGS="-O2"
    CXXFLAGS=""
  elif test "x$with_opt" = x0; then
    AC_MSG_RESULT([${with_opt}])
    OPTFLAGS="-O0"
    CXXFLAGS=""
  elif test "x$with_opt" = x1; then
    AC_MSG_RESULT([${with_opt}])
    OPTFLAGS="-O1"
    CXXFLAGS=""
  elif test "x$with_opt" = x2; then
    dnl Default case is explicitly specified
    AC_MSG_RESULT([${with_opt}])
    OPTFLAGS="-O2"
    CXXFLAGS=""
  elif test "x$with_opt" = x3; then
    AC_MSG_RESULT([${with_opt}])
    OPTFLAGS="-O3"
    CXXFLAGS=""
  else
    AC_MSG_RESULT([${with_opt}])
    AC_MSG_NOTICE([Unrecognized optimization level, using default '-02'.])
    dnl OPTFLAGS="-O2"
    dnl CXXFLAGS=""
  fi
  AC_SUBST(OPTFLAGS)
  AM_CONDITIONAL([WITH_CUSTOM_OPT_LEVEL],[test "x$OPTFLAGS" != "x"])

])
