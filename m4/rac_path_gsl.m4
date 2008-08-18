dnl @synopsis RAC_PATH_GSL
dnl
dnl Allow path to GSL to be specified on the configure line.
dnl

AC_DEFUN([RAC_PATH_GSL],[
  AC_ARG_WITH([gsl],
              AC_HELP_STRING([--with-gsl=DIR], [use GSL (default is YES),
                specify the root directory for GSL library (optional)]),
                [ if test "$withval" = "no"; then
                    want_gsl="no"
                  elif test "$withval" = "yes"; then
                    want_gsl="yes"
                    ac_gsl_path=""
                  else
                    want_gsl="yes"
                    ac_gsl_path="$withval"
                  fi
                ],[want_gsl="yes"])

  AC_MSG_CHECKING([whether with third-party library, GSL, is wanted])
  AC_MSG_RESULT([${want_gsl}])

  if test "x$want_gsl" = xyes; then
    AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])
    AC_MSG_CHECKING([for GSL directory])

    if test "$ac_gsl_path" != ""; then
      AC_MSG_RESULT([${ac_gsl_path}])
      if test -d "${ac_gsl_path}"; then
        GSL_CPPFLAGS="-I$ac_gsl_path"
        GSL_LDFLAGS="-L$ac_gsl_path"

        AC_MSG_NOTICE(GSL_CPPFLAGS: $GSL_CPPFLAGS)
        AC_SUBST(GSL_CPPFLAGS)
        AC_SUBST(GSL_LDFLAGS)
      else
        AC_MSG_ERROR([the specified GSL basedir, ${ac_gsl_path} does not exist])
      fi
    else
      dnl FLAGS for default case (SVN tree) are easily set in the Makefile.am
      AC_MSG_RESULT([using default GSL directory in the local SVN tree])
      AC_CONFIG_SUBDIRS([packages/gsl])
    fi
  fi
  AM_CONDITIONAL([WITH_GSL],[test "x$with_gsl" != xno])
])
