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
    AC_DEFINE([HAVE_DFFTPACK],[1], [Macro to handle code which depends on DFFTPACK.])
    DFFTPACK_LDFLAGS="-L`pwd`/packages/dfftpack"
    AC_SUBST(DFFTPACK_LDFLAGS)
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
  if test "x$with_lhs" = xyes -a "x$enable_f90" = xyes; then
    AC_DEFINE([HAVE_LHS],[1],[Macro to handle code which depends on LHS.])
    AC_CONFIG_SUBDIRS([packages/LHS])
    LHS_LDFLAGS="-L`pwd`/packages/LHS"
    AC_SUBST(LHS_LDFLAGS)
  fi
  AM_CONDITIONAL([WITH_LHS],
                 [test "x$with_lhs" = xyes])

  dnl Teuchos package checks.
  AC_ARG_WITH([teuchos],
              AC_HELP_STRING([--with-teuchos=<dir>],
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
              AC_HELP_STRING([--with-gsl=<dir>],
                             [use GSL (default is yes), specify the root
                              directory for GSL library (optional)]),
              [],[with_gsl="yes"])

  acx_local_gsl=no
  case $with_gsl in
  no)
    AC_MSG_NOTICE([NOT building with GSL!])
    GSL_ROOT=""
    ;;

  system)
    AC_MSG_CHECKING([for system GSL])
    AC_CHECK_LIB([gslcblas],[cblas_dgemm])
    AC_CHECK_LIB([gsl],[gsl_ran_gamma_pdf])
    AC_CHECK_HEADERS([gsl/gsl_version.h])
    AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])
    AC_MSG_NOTICE([Using system GSL!])
    GSL_ROOT=""
    ;;

  dnl For yes, check GSL_ROOT, otherwise fallback to GSL from SVN checkout
  yes | "")
    AC_MSG_CHECKING([for GSL])
    if test -n "$GSL_ROOT" -a -d "$GSL_ROOT"; then

      AC_MSG_NOTICE([GSL_ROOT is set by the caller])
      dnl save_CPPFLAGS=$CPPFLAGS
      dnl CPPFLAGS="-I$GSL_ROOT"
      dnl AC_CHECK_HEADERS([gsl/gsl_version.h],acx_external_gsl=yes,,"-I$GSL_ROOT")
      dnl CPPFLAGS=$save_CPPFLAGS
      AC_MSG_RESULT([using GSL in GSL_ROOT: $GSL_ROOT])
      AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])

    elif test -d `pwd`/packages/gsl; then

      dnl use SVN checkout of GSL and instruct subpackages to do so as well
      GSL_ROOT=`pwd`/packages/gsl
      acx_local_gsl=yes
      AC_CONFIG_SUBDIRS([packages/gsl])

      dnl no tests to perform since GSL has yet to be built; trust the SVN co
      AC_MSG_RESULT([using GSL in $GSL_ROOT])
      AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])

    else
      AC_MSG_RESULT([could not find GSL directory.])
    fi
    ;;

  dnl Otherwise, user should have provided an explicit path to GSL
  *)
    AC_MSG_CHECKING([for GSL])
    GSL_ROOT=$withval
    if test -n "$GSL_ROOT" -a -d "$GSL_ROOT"; then

      AC_MSG_RESULT([using: $GSL_ROOT])

      dnl For now, trust the user (future enhancement would provide addl checks)
      dnl AC_CHECK_HEADERS([gsl/gsl_version.h])
      AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])
      AC_MSG_NOTICE([Using the GSL specified on the configure line!])

    else
      AC_MSG_ERROR([could not locate $GSL_ROOT])
    fi
    ;;
  esac

  dnl Do not export GSL build variables in certain cases (i.e. no or system)
  if test -n "$GSL_ROOT"; then
    GSL_CPPFLAGS="-I$GSL_ROOT"
    GSL_LDFLAGS="-L$GSL_ROOT"

    AC_SUBST(GSL_ROOT)
    AC_SUBST(GSL_CPPFLAGS)
    AC_SUBST(GSL_LDFLAGS)
  fi

  dnl WITH_GSL controls CPPFLAGS in src/Makefile.am
  AM_CONDITIONAL([WITH_GSL], [test "x$with_gsl" = xsystem -o -n "$GSL_ROOT" ])
  dnl BUILD_GSL controls build of Pecos local gsl in packages/Makefile.am
  AM_CONDITIONAL([BUILD_GSL], [test "x$acx_local_gsl" = xyes])

])
