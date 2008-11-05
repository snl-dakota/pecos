## Process this file with automake to produce 'Makefile.in'.

AM_CPPFLAGS  = -I${top_srcdir}/src @TEUCHOS_CPPFLAGS@
TPL_LDFLAGS  = @TEUCHOS_LDFLAGS@
TPL_LIBS     = -lteuchos

if WITH_GSL
AM_CPPFLAGS += @GSL_CPPFLAGS@
TPL_LDFLAGS += @GSL_LDFLAGS@
TPL_LIBS    += -lgsl
endif

if WITH_LHS
TPL_LDFLAGS += @LHS_LDFLAGS@
TPL_LIBS    += -llhs
endif

if WITH_DFFTPACK
TPL_LDFLAGS += @DFFTPACK_LDFLAGS@
TPL_LIBS    += -ldfftpack
endif


# Tests to run when 'gmake check' is run
TESTS = test_pecos_ifft_g test_pecos_ifft_sd

# Programs to link when 'gmake' is run
noinst_PROGRAMS = pecos_ifft_g pecos_ifft_sd


# Complete these for each specific test.
pecos_ifft_g_SOURCES = pecos_ifft_g.cpp
pecos_ifft_g_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_g_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

# Complete these for each specific test.
pecos_ifft_sd_SOURCES = pecos_ifft_sd.cpp
pecos_ifft_sd_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_sd_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

# Define this variable if you want to install test header files.
# subpkginclude_HEADERS =
