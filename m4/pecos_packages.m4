dnl Packages

AC_DEFUN([PECOS_PACKAGES],[
  dnl FFT packages - BOTH dfftpack and fftw will be built
  AC_ARG_WITH([fft],AS_HELP_STRING([--without-fft],
              [turn FFT support off]),[with_fft=$withval],[with_fft=no])
  if test "x$with_fft" = xyes; then
    dnl AC_CONFIG_SUBDIRS([packages/dfftpack])
    dnl AC_CONFIG_SUBDIRS([packages/fftw])
    AC_DEFINE([HAVE_FFT],[1], [Macro to handle code which depends on FFT.])
  fi
  AM_CONDITIONAL([WITH_FFT],[test "x$with_fft" = xyes])

  dnl DFFTPACK package checks.
  AC_ARG_WITH([dfftpack],AS_HELP_STRING([--without-dfftpack],
              [turn DFFTPACK support off]),[with_dfftpack=$withval],[with_dfftpack=yes])
  if test "x$with_dfftpack" = xyes; then
    AC_CONFIG_SUBDIRS([packages/dfftpack])
    AC_DEFINE([PECOS_DFFTPACK],[1], [Macro to handle code which depends on DFFTPACK.])
  fi
  AM_CONDITIONAL([WITH_DFFTPACK],[test "x$with_dfftpack" = xyes])

  dnl FFTW package checks.
  AC_ARG_WITH([fftw],AS_HELP_STRING([--without-fftw],
              [turn FFTW support off]),[with_fftw=$withval],[with_fftw=no])
  if test "x$with_fftw" = xyes; then
    AC_CONFIG_SUBDIRS([packages/fftw])
    AC_DEFINE([HAVE_FFTW],[1], [Macro to handle code which depends on FFTW.])
  fi
  AM_CONDITIONAL([WITH_FFTW],[test "x$with_fftw" = xyes])

  dnl LHS package checks.
  AC_ARG_WITH([lhs],AS_HELP_STRING([--without-lhs],[turn LHS support off]),
	      [with_lhs=$withval],[with_lhs=yes])
  if test "x$with_lhs" = xyes; then
    AC_DEFINE([PECOS_LHS],[1],[Macro to handle code which depends on LHS.])
    #if test "x$enable_f90" = xyes; then
      AC_CONFIG_SUBDIRS([packages/LHS])
    #fi
  fi
  AM_CONDITIONAL([WITH_LHS],
                 [test "x$with_lhs" = xyes])

  dnl Teuchos package checks.
  AC_ARG_WITH([teuchos],
              AC_HELP_STRING([--with-teuchos=DIR],
                             [use Teuchos (default is yes), specify the root
                              directory for Teuchos library]),
              [],[with_teuchos="yes"])
  acx_external_teuchos=no
  case $with_teuchos in
  dnl Pecos depends on Teuchos UNCONDITIONALLY
  no)
    AC_MSG_ERROR([Pecos cannot be configured without Teuchos. Please specify
                 --with-teuchos OR provide a path to a prebuilt Teuchos.])
    ;;

  dnl For yes, check TEUCHOS_ROOT, otherwise fallback to local Teuchos
  yes | "")
    AC_MSG_CHECKING([for Teuchos])
    if test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT"; then

      acx_external_teuchos=yes
      AC_MSG_RESULT([using Teuchos in TEUCHOS_ROOT: $TEUCHOS_ROOT])

    elif test -d `pwd`/packages/teuchos; then

      dnl use local teuchos and instruct subpackages to do so as well
      export TEUCHOS_ROOT=`pwd`/packages/teuchos
      acx_external_teuchos=no
      AC_CONFIG_SUBDIRS([packages/teuchos])
      AC_MSG_RESULT([using local Teuchos in $TEUCHOS_ROOT])

    else
      AC_MSG_RESULT([could not find Teuchos directory.])
    fi
    ;;

  dnl Otherwise, user should have provided an explicit path to Teuchos
  *)
    AC_MSG_CHECKING([for Teuchos])
    TEUCHOS_ROOT=$withval
    if test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT"; then

      acx_external_teuchos=yes
      AC_MSG_RESULT([using: $TEUCHOS_ROOT])

    else
      AC_MSG_ERROR([could not locate $TEUCHOS_ROOT])
    fi
    ;;
  esac

  TEUCHOS_CPPFLAGS="-I$TEUCHOS_ROOT/src"
  TEUCHOS_LDFLAGS="-L$TEUCHOS_ROOT/src"

  AC_SUBST(TEUCHOS_ROOT)
  AC_SUBST(TEUCHOS_CPPFLAGS)
  AC_SUBST(TEUCHOS_LDFLAGS)

  AM_CONDITIONAL([WITH_ALT_EXTERNAL_TEUCHOS],
                 [test "x$acx_external_teuchos" = xyes])

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
