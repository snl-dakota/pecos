dnl Packages

AC_DEFUN([PECOS_PACKAGES],[
  dnl BOOST package package check.
  AC_ARG_WITH([boost],AS_HELP_STRING([--without-boost],
              [turn Boost support off]),[with_boost=$withval],[with_boost=yes])
  if test "x$with_boost" = xyes; then
    dnl AC_CONFIG_SUBDIRS([packages/boost])
    dnl AC_DEFINE([HAVE_BOOST],[1],[Macro to handle code which depends on Boost.])
    BOOST_CPPFLAGS="-I`pwd`/packages"
    AC_SUBST(BOOST_CPPFLAGS)
  else
    AC_MSG_NOTICE([could not find boost directory!])
    AC_MSG_ERROR([please ensure your Pecos distribution includes boost header
                  files in: <pecos_root>/packages/boost])
  fi
  AM_CONDITIONAL([WITH_BOOST],[test "x$with_boost" = xyes])

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
              [turn FFTW support off]),[with_fftw=$withval],[with_fftw=yes])
  if test "x$with_fftw" = xyes; then
    AC_CONFIG_SUBDIRS([packages/fftw])
    AC_DEFINE([HAVE_FFTW],[1], [Macro to handle code which depends on FFTW.])
    FFTW_CPPFLAGS="-I`pwd`/packages/fftw/api"
    FFTW_LDFLAGS="-L`pwd`/packages/fftw"
    AC_SUBST(FFTW_CPPFLAGS)
    AC_SUBST(FFTW_LDFLAGS)
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
  AM_CONDITIONAL([WITH_LHS],[test "x$with_lhs" = xyes -a "x$enable_f90" = xyes])

  dnl Teuchos package checks.
  AC_ARG_WITH([teuchos],
              AC_HELP_STRING([--with-teuchos=<dir>],
                             [use Teuchos (default is yes), specify the root
                              directory for Teuchos library]),
              [],[with_teuchos="yes"])
  acx_local_teuchos=no
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

      AC_MSG_RESULT([using Teuchos in TEUCHOS_ROOT: $TEUCHOS_ROOT])

    elif test -d `pwd`/packages/teuchos; then

      dnl use local teuchos and instruct subpackages to do so as well
      export TEUCHOS_ROOT=`pwd`/packages/teuchos
      acx_local_teuchos=yes
      AC_CONFIG_SUBDIRS([packages/teuchos])
      AC_MSG_RESULT([using local Teuchos in $TEUCHOS_ROOT])

    else
      AC_MSG_NOTICE([could not find Teuchos directory.])
      AC_MSG_NOTICE([need help locating teuchos!])
      AC_MSG_ERROR([PLEASE PROVIDE full path to teuchos, --with-teuchos=<DIR>])
    fi
    ;;

  dnl Otherwise, user should have provided an explicit path to Teuchos
  *)
    AC_MSG_CHECKING([for Teuchos])
    TEUCHOS_ROOT=$withval
    if test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT"; then
      AC_MSG_RESULT([using: $TEUCHOS_ROOT])
    else
      AC_MSG_ERROR([could not locate $TEUCHOS_ROOT])
    fi
    ;;
  esac

  dnl Finally, check for INSTALLED Teuchos vs. BUILT, but NOT-installed Teuchos
  if test -d "$TEUCHOS_ROOT/include" -a -d "$TEUCHOS_ROOT/lib"; then
    AC_MSG_NOTICE([Found an INSTALLED teuchos!])
    TEUCHOS_CPPFLAGS="-I$TEUCHOS_ROOT/include"
    TEUCHOS_LDFLAGS="-L$TEUCHOS_ROOT/lib"
  elif test -d "$TEUCHOS_ROOT/src"; then
    TEUCHOS_CPPFLAGS="-I$TEUCHOS_ROOT/src"
    TEUCHOS_LDFLAGS="-L$TEUCHOS_ROOT/src"
  else
    AC_MSG_ERROR([could not find Teuchos library relative to $TEUCHOS_ROOT.])
  fi

  AC_SUBST(TEUCHOS_ROOT)
  AC_SUBST(TEUCHOS_CPPFLAGS)
  AC_SUBST(TEUCHOS_LDFLAGS)

  AM_CONDITIONAL([BUILD_TEUCHOS], [test "x$acx_local_teuchos" = xyes])
])
