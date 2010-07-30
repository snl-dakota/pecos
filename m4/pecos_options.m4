dnl PECOS Options

AC_DEFUN([PECOS_OPTIONS],[
  dnl Debug option check.
  AC_ARG_ENABLE([debug],
                AS_HELP_STRING([--enable-debug],[turn debug support on]),
                [enable_debug=$enableval],[enable_debug=no])
  dnl AC_MSG_CHECKING([debug])
  if test "x$enable_debug" = xyes; then
    dnl AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_PECOS_DEBUG],[1],[Macro to enable debug in Pecos.])
    DEBUGFLAGS="-g"
  else
    dnl AC_MSG_RESULT([no])
    AC_DEFINE([BOOST_DISABLE_ASSERTS],[1],[Macro to disable Boost asserts.])
    DEBUGFLAGS=""
  fi
  AC_SUBST(DEBUGFLAGS)
  AM_CONDITIONAL([WITH_DEBUG_ENABLED],[test "x$enable_debug" = xyes])

  dnl Tests option check.
  AC_ARG_ENABLE([tests],
                AS_HELP_STRING([--disable-tests],[do not build unit tests]),
                [enable_tests=$enableval],[enable_tests=yes])
  if test "x$enable_debug" = xyes; then
    AC_CONFIG_SUBDIRS([test])
  fi
  AM_CONDITIONAL([ENABLE_TESTS],[test "x$enable_tests" = xyes])
])
