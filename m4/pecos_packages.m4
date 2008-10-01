dnl Packages

AC_DEFUN([PECOS_PACKAGES],[
  dnl FFT packages - BOTH dfftpack and fftw will be built
  AC_ARG_WITH([fft],AS_HELP_STRING([--without-fft],
              [turn FFT support off]),[with_fft=$withval],[with_fft=no])
  if test "x$with_fft" = xyes; then
    AC_CONFIG_SUBDIRS([packages/dfftpack])
    dnl AC_CONFIG_SUBDIRS([packages/fftw])
    dnl AC_DEFINE([HAVE_FFT],[1], [Macro to handle code which depends on FFT.])
  fi
  AM_CONDITIONAL([WITH_FFT],[test "x$with_fft" = xyes])

  dnl DFFTPACK package checks.
  AC_ARG_WITH([dfftpack],AS_HELP_STRING([--without-dfftpack],
              [turn DFFTPACK support off]),[with_dfftpack=$withval],[with_dfftpack=yes])
  if test "x$with_dfftpack" = xyes; then
    dnl AC_CONFIG_SUBDIRS([packages/dfftpack])
    AC_DEFINE([HAVE_DFFTPACK],[1], [Macro to handle code which depends on DFFTPACK.])
  fi
  AM_CONDITIONAL([WITH_DFFTPACK],[test "x$with_dfftpack" = xyes])

  dnl FFTW package checks.
  AC_ARG_WITH([fftw],AS_HELP_STRING([--without-fftw],
              [turn FFTW support off]),[with_fftw=$withval],[with_fftw=no])
  if test "x$with_fftw" = xyes; then
    dnl AC_CONFIG_SUBDIRS([packages/fftw])
    AC_DEFINE([HAVE_FFTW],[1], [Macro to handle code which depends on FFTW.])
  fi
  AM_CONDITIONAL([WITH_FFTW],[test "x$with_fftw" = xyes])

  dnl Teuchos package checks.
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
    dnl WJB TRASH AC_MSG_RESULT([${withval}])

  dnl Teuchos dependency can be managed with an alternate, external path setting
  AM_CONDITIONAL([WITH_ALT_EXTERNAL_TEUCHOS],[test "x$want_teuchos" = xyes -a \
                                                   -d "${ac_teuchos_path}"])

  if test "x$want_teuchos" = xno; then
    dnl Pecos depends on Teuchos UNCONDITIONALLY
    AC_MSG_ERROR([Pecos cannot be configured without Teuchos. Please specify --with-teuchos=yes OR provide a DIR path to a prebuilt Teuchos])
  else
    AC_MSG_CHECKING([for Teuchos directory])

    if test "$ac_teuchos_path" = ""; then
      dnl OVERRIDE 'yes' case with directory path resulting from SVN checkout

      ac_teuchos_path=`cd ${ac_top_builddir}. && pwd`/packages/teuchos
      dnl AC_MSG_NOTICE(ac_teuchos_path: $ac_teuchos_path)

      AC_MSG_RESULT([the default Teuchos basedir, ${ac_teuchos_path}, will be configured and used for this Pecos build])
      AC_CONFIG_SUBDIRS([packages/teuchos])
    else
      AC_MSG_RESULT([${ac_teuchos_path}])
    fi

    if test -d "${ac_teuchos_path}"; then
      TEUCHOS_CPPFLAGS="-I$ac_teuchos_path/src"
      TEUCHOS_LDFLAGS="-L$ac_teuchos_path/src"

      AC_MSG_NOTICE(TEUCHOS_CPPFLAGS: $TEUCHOS_CPPFLAGS)
      AC_SUBST(TEUCHOS_CPPFLAGS)
      AC_SUBST(TEUCHOS_LDFLAGS)
    else
      AC_MSG_ERROR([the specified Teuchos basedir, ${ac_teuchos_path} does not exist])
    fi
  fi

  dnl GSL package checks.
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
