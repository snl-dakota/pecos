dnl Packages

AC_DEFUN([PECOS_PACKAGES],[
  AC_ARG_WITH([boost],
              AC_HELP_STRING([--with-boost=DIR],
                             [use Boost headers in specified DIR]),
              [],[with_boost="yes"])
  case $with_boost in
  no) 
    dnl BOOST package is needed unconditionally.
    dnl Pecos provides a header-only subset of the Boost 1.40 release.
    AC_MSG_ERROR([PECOS cannot be configured without Boost. Please specify
                  directory to boost headers OR simply, --with-boost=yes
                  to get the default path to the PECOS provided boost subset.])
    ;;      

  dnl For yes, check BOOST_ROOT, otherwise fallback to local Boost
  yes | "")
    AC_MSG_CHECKING([for Boost])
    if test -n "$BOOST_ROOT" -a -d "$BOOST_ROOT"; then

      AC_MSG_RESULT([using Boost in BOOST_ROOT: $BOOST_ROOT])

    elif test -d `pwd`/packages/boost; then

      dnl Use local Boost and instruct subpackages to do so as well
      export BOOST_ROOT=`pwd`/packages/boost

      dnl Nothing to config/build; Pecos provides a header-only subset
      dnl AC_CONFIG_SUBDIRS([packages/boost])
      AC_MSG_RESULT([using local Boost in $BOOST_ROOT])

    else
      AC_MSG_NOTICE([could not find Boost directory.])
      AC_MSG_NOTICE([need help locating boost!])
      AC_MSG_ERROR([PLEASE PROVIDE full path to boost, --with-boost=<DIR>])
    fi
    ;;

  dnl Otherwise, user should have provided an explicit path to Boost
  *)
    AC_MSG_CHECKING([for specified Boost])
    BOOST_ROOT=$withval
    if test -n "$BOOST_ROOT" -a -d "$BOOST_ROOT"; then
      AC_MSG_RESULT([using: $BOOST_ROOT])
    else
      AC_MSG_ERROR([could not locate $BOOST_ROOT])
    fi
    ;;

  esac

  boost_version_req=103700
  boost_version=`grep 'define BOOST_VERSION 1' $BOOST_ROOT/boost/version.hpp | cut -d' ' -f3`
  if test $boost_version -ge $boost_version_req; then
    AC_MSG_RESULT([Boost meets PECOS min version: $boost_version >= $boost_version_req])
  else
    AC_MSG_ERROR([OLD Boost: $boost_version needs to be $boost_version_req or later])
  fi

  BOOST_CPPFLAGS="-I$BOOST_ROOT"
  AC_SUBST(BOOST_CPPFLAGS)
  AC_ARG_VAR(BOOST_ROOT, [Path to header-only subset of Boost, a C++ foundation package])


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
              [turn support for GPL package, FFTW off]),[with_fftw=$withval],[with_fftw=no])
  if test "x$with_fftw" = xyes; then
    AC_CONFIG_SUBDIRS([packages/fftw])
    AC_DEFINE([HAVE_FFTW],[1], [Macro to handle code which depends on FFTW.])
    FFTW_CPPFLAGS="-I`pwd`/packages/fftw/api"
    FFTW_LDFLAGS="-L`pwd`/packages/fftw"
    AC_SUBST(FFTW_CPPFLAGS)
    AC_SUBST(FFTW_LDFLAGS)
    AC_MSG_NOTICE([NOTE: your PECOS build includes the binary, GPL library, FFTW!])
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

  dnl VPISparseGrid package checks.
  AC_ARG_WITH([sparsegrid],AS_HELP_STRING([--without-sparsegrid],
	      [turn VPISparseGrid support off]),[with_sparsegrid=$withval],
	      [with_sparsegrid=yes])
  if test "x$with_sparsegrid" = xyes; then
    AC_DEFINE([HAVE_SPARSE_GRID],[1],
	      [Macro to handle code which depends on VPISparseGrid.])
    AC_CONFIG_SUBDIRS([packages/VPISparseGrid])
    dnl The following macros are used in Pecos/test/Makefile.am
    SPARSEGRID_CPPFLAGS="-I`pwd`/packages/VPISparseGrid/src"
    SPARSEGRID_LDFLAGS="-L`pwd`/packages/VPISparseGrid/src"
    AC_SUBST(SPARSEGRID_CPPFLAGS)
    AC_SUBST(SPARSEGRID_LDFLAGS)
  fi
  AM_CONDITIONAL([WITH_SPARSE_GRID],[test "x$with_sparsegrid" = xyes])


  dnl -------------------------
  dnl Teuchos include DIR check
  dnl -------------------------
  AC_ARG_WITH([teuchos-include],
              AC_HELP_STRING([--with-teuchos-include=DIR],
                             [use Teuchos headers in specified include DIR]),
              [],[with_teuchos_include="yes"])
  case $with_teuchos_include in
  
  no) 
    AC_MSG_ERROR([Usage: --with-teuchos-include=DIR; without not allowed])
    ;;      

  yes | "")
    AC_MSG_CHECKING([for Teuchos include dir via environment variable])
    if test -n "$TEUCHOS_INCLUDE" -a -d "$TEUCHOS_INCLUDE"; then
      dnl could check for existence of header
      AC_MSG_RESULT([using Teuchos headers in TEUCHOS_INCLUDE: $TEUCHOS_INCLUDE])   
    else
      AC_MSG_NOTICE([Teuchos include dir not specfied; will find later.])
    fi
    ;;

  *)
    AC_MSG_CHECKING([for specified Teuchos include DIR])
    TEUCHOS_INCLUDE=$withval
    if test -n "$TEUCHOS_INCLUDE" -a -d "$TEUCHOS_INCLUDE"; then
      export TEUCHOS_INCLUDE=$withval
      AC_MSG_RESULT([using Teuchos headers in $TEUCHOS_INCLUDE])
    else
      AC_MSG_ERROR([could not locate specified Teuchos include dir $TEUCHOS_INCLUDE])
    fi
    ;;

  esac

  AC_ARG_VAR(TEUCHOS_INCLUDE,
             [Path to headers for Teuchos, OO Numerics foundation library])


  dnl -------------------------
  dnl Teuchos lib DIR check
  dnl -------------------------
  AC_ARG_WITH([teuchos-lib],
              AC_HELP_STRING([--with-teuchos-lib=DIR],
                             [use Teuchos libraries in specified lib DIR]),
              [],[with_teuchos_lib="yes"])

  case $with_teuchos_lib in

  no)
    AC_MSG_ERROR([Usage: --with-teuchos-lib=DIR; without not allowed])
    ;;

  yes | "")

    AC_MSG_CHECKING([for Teuchos lib dir via environment variable])
    if test -n "$TEUCHOS_LIB" -a -d "$TEUCHOS_LIB"; then
      dnl cannot check for existence of lib in this case,
      dnl since might not be installed
      AC_MSG_RESULT([using Teuchos libs in TEUCHOS_LIB: $TEUCHOS_LIB])
    else
      AC_MSG_NOTICE([Teuchos lib dir not specfied; will find later.])
    fi
    ;;

  *)
    AC_MSG_CHECKING([for specified Teuchos lib DIR])
    TEUCHOS_LIB=$withval
    if test -n "$TEUCHOS_LIB" -a -d "$TEUCHOS_LIB"; then
      export TEUCHOS_LIB=$withval
      AC_MSG_RESULT([using Teuchos libs in $TEUCHOS_LIB])
    else
      AC_MSG_ERROR([could not locate specified Teuchos lib dir $TEUCHOS_LIB])
    fi
    ;;

  esac

  AC_ARG_VAR(TEUCHOS_LIB,
             [Path to libraries for Teuchos OO Numerics foundation library])

  dnl Teuchos package checks: this will perform the final checks to
  dnl resolve all Teuchos-related configure options and make file
  dnl include/lib settings
  AC_ARG_WITH([teuchos],
              AC_HELP_STRING([--with-teuchos=DIR],
                             [use Teuchos (default is yes), optionally
                              specify the root Teuchos directory containing
                              src or include/lib]),
              [],[with_teuchos="yes"])

  dnl should also include test for directory either here or above
  acx_valid_inc=0
  if test -n "$TEUCHOS_INCLUDE" -a -d "$TEUCHOS_INCLUDE"; then
    acx_valid_inc=1
  fi
  acx_valid_lib=0
  if test -n "$TEUCHOS_LIB" -a -d "$TEUCHOS_LIB"; then
    acx_valid_lib=1
  fi

  dnl this is cludgy, but not sure how to do booleans nor xor
  dnl is shell arithmetic safe here?
  if test $acx_valid_inc -eq 1 -a $acx_valid_lib -eq 0; then
    AC_MSG_ERROR([Must specify both or neither of Teuchos include/lib, not exactly one.])
  fi
  if test $acx_valid_inc -eq 0 -a $acx_valid_lib -eq 1; then
    AC_MSG_ERROR([Must specify both or neither of Teuchos include/lib, not exactly one.])
  fi

  acx_local_teuchos=no
  case $with_teuchos in
  dnl DAKOTA depends on Teuchos UNCONDITIONALLY
  no)
    AC_MSG_ERROR([DAKOTA cannot be configured without Teuchos. Please use a
                 --with-teuchos option OR provide path to a prebuilt Teuchos.])
    ;;

  dnl For yes, check environment variables, otherwise fallback to local Teuchos
  yes | "")
    AC_MSG_CHECKING([for Teuchos])
    if test $acx_valid_inc -eq 1 -a $acx_valid_lib -eq 1; then
      AC_MSG_RESULT([using previously specified Teuchos include/lib])
    elif test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT"; then
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
    if test $acx_valid_inc -eq 1 -a $acx_valid_lib -eq 1; then
      AC_MSG_ERROR([Teuchos include/lib and root DIR are mutually exclusive])
    fi
    AC_MSG_CHECKING([for specified Teuchos])
    TEUCHOS_ROOT=$withval
    if test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT"; then
      AC_MSG_RESULT([using: $TEUCHOS_ROOT])
    else
      AC_MSG_ERROR([could not locate $TEUCHOS_ROOT])
    fi
    ;;

  esac

  if test $acx_valid_inc -eq 1 -a $acx_valid_lib -eq 1; then

    TEUCHOS_CPPFLAGS="-I$TEUCHOS_INCLUDE -I$TEUCHOS_LIB"
    TEUCHOS_LDFLAGS="-L$TEUCHOS_LIB"

  dnl Check for INSTALLED Teuchos vs. BUILT, but NOT-installed Teuchos
  elif test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT/include" -a -d "$TEUCHOS_ROOT/lib"; then

    AC_MSG_NOTICE([Found an INSTALLED teuchos!])
    TEUCHOS_CPPFLAGS="-I$TEUCHOS_ROOT/include"
    TEUCHOS_LDFLAGS="-L$TEUCHOS_ROOT/lib"

  elif test -n "$TEUCHOS_ROOT" -a -d "$TEUCHOS_ROOT/src"; then

    AC_MSG_NOTICE([Found a source teuchos!])
    TEUCHOS_CPPFLAGS="-I$TEUCHOS_ROOT/src"
    TEUCHOS_LDFLAGS="-L$TEUCHOS_ROOT/src"

  else

    AC_MSG_ERROR([could not find Teuchos library relative to root nor include/lib.])

  fi

  AC_MSG_NOTICE([Final Teuchos config CPPFLAGS: $TEUCHOS_CPPFLAGS LDFLAGS: $TEUCHOS_LDFLAGS])

  AC_ARG_VAR(TEUCHOS_ROOT, [Path to Teuchos, OO Numerics foundation library])
  AC_SUBST(TEUCHOS_CPPFLAGS)
  AC_SUBST(TEUCHOS_LDFLAGS)

  AM_CONDITIONAL([BUILD_TEUCHOS], [test "x$acx_local_teuchos" = xyes])
])
