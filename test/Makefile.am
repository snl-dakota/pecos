# WJB - ToDo: look into the possibility of "sharing" across all Pecos Makefiles
####include $(top_srcdir)/Pecos.am.mk

# This definition is needed to make full use of Pecos.am.mk, e.g init AM_CPPFLAGS (set to -I${top_srcdir}/src so += operator to concat CPPFLAGS)
####THISDIR = Pecos/test

AM_CPPFLAGS  = -I${top_srcdir}/src @TEUCHOS_CPPFLAGS@
AM_LDFLAGS   = @TEUCHOS_LDFLAGS@
TPL_LIBS     = -lteuchos

if WITH_GSL
AM_CPPFLAGS += @GSL_CPPFLAGS@
AM_LDFLAGS  += @GSL_LDFLAGS@
TPL_LIBS    += -lgsl
endif

if WITH_LHS
AM_LDFLAGS  += @LHS_LDFLAGS@
TPL_LIBS    += -llhs
endif


# Tests to run when 'gmake check' is run
TESTS = test_pecos_ifft_g test_pecos_ifft_sd

# Programs to link when 'gmake' is run
noinst_PROGRAMS = pecos_ifft_g pecos_ifft_sd


# Complete these for each specific test.
pecos_ifft_g_SOURCES = pecos_ifft_g.cpp
pecos_ifft_g_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_g_LDADD = -L${top_builddir} -lpecos \
    -L${top_builddir}/packages/dfftpack -ldfftpack $(TPL_LIBS) \
    @LAPACK_LIBS@ @BLAS_LIBS@ @FCLIBS@

# Complete these for each specific test.
pecos_ifft_sd_SOURCES = pecos_ifft_g.cpp
pecos_ifft_sd_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_sd_LDADD = -L${top_builddir} -lpecos \
    -L${top_builddir}/packages/dfftpack -ldfftpack $(TPL_LIBS) \
    @LAPACK_LIBS@ @BLAS_LIBS@ @FCLIBS@

# Define this variable if you want to install test header files.
# subpkginclude_HEADERS =
