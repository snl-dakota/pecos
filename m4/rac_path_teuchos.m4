dnl @synopsis RAC_PATH_TEUCHOS
dnl
dnl Allow path to Teuchos to be specified on the configure line.
dnl

AC_DEFUN([RAC_PATH_TEUCHOS],[
  AC_ARG_WITH([teuchos],
              AC_HELP_STRING([--with-teuchos=DIR], [use Teuchos (default is YES),
                specify the root directory for Teuchos library (RECOMMENDED)]),
                [ if test "$withval" = "no"; then
                    want_teuchos="no"
                  elif test "$withval" = "yes"; then
                    want_teuchos="yes"
                    ac_teuchos_path=""
                  else
                    want_teuchos="yes"
                    ac_teuchos_path="$withval"
                  fi
                ],[want_teuchos="yes"])

  AC_MSG_CHECKING([whether with third-party library, Teuchos, is wanted])
  AC_MSG_RESULT([${want_teuchos}])

  if test "x$want_teuchos" = xyes; then
    dnl Pecos depends on Teuchos UNCONDITIONALLY
    dnl AC_DEFINE([TPL_TEUCHOS],[1],[Macro to handle code which depends on Teuchos.])
    AC_MSG_CHECKING([for Teuchos directory])

    if test "$ac_teuchos_path" != ""; then
      AC_MSG_RESULT([${ac_teuchos_path}])
      if test -d "${ac_teuchos_path}"; then
        TEUCHOS_CPPFLAGS="-I$ac_teuchos_path/src"
        TEUCHOS_LDFLAGS="-L$ac_teuchos_path/src"

        AC_MSG_NOTICE(TEUCHOS_CPPFLAGS: $TEUCHOS_CPPFLAGS)
        AC_SUBST(TEUCHOS_CPPFLAGS)
        AC_SUBST(TEUCHOS_LDFLAGS)
      else
        AC_MSG_ERROR([the specified Teuchos basedir, ${ac_teuchos_path} does not exist])
      fi
    else
      dnl FLAGS for default case (SVN tree) are easily set in the Makefile.am
      AC_MSG_RESULT([using default Teuchos directory in the local SVN tree])
      AC_CONFIG_SUBDIRS([packages/teuchos])
    fi

  else
    AC_MSG_ERROR([Pecos cannot be configured without Teuchos. Please specify --with-teuchos=yes OR provide a DIR path to Teuchos])
  fi

  AM_CONDITIONAL([WITH_ALT_EXTERNAL_TEUCHOS],[test "x$want_teuchos" = xyes && test -d "${ac_teuchos_path}"])
])
