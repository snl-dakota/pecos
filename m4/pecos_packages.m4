dnl Packages

AC_DEFUN([PECOS_PACKAGES],[
  
  dnl -------------------
  dnl Boost package check
  dnl -------------------
  AX_BOOST_BASE([1.37], [required], [$srcdir/packages/boost], [packages/boost])

  if test "$ac_boost_build_tpl" = "yes"; then
    AC_MSG_NOTICE([will build bundled boost TPL.])
    AC_CONFIG_SUBDIRS([packages/boost])
  else
    AC_MSG_NOTICE([skipping bundled boost TPL.])
  fi

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
  fi
  AM_CONDITIONAL([WITH_SPARSE_GRID],[test "x$with_sparsegrid" = xyes])

  dnl ---------------------
  dnl Teuchos package check
  dnl ---------------------
  PECOS_AC_TEUCHOS()

])
