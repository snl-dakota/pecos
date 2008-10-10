# WJB - ToDo: look into the possibility of "sharing" across all Pecos Makefiles
####include $(top_srcdir)/Pecos.am.mk

# This definition is needed to make full use of Pecos.am.mk, e.g init AM_CPPFLAGS (set to -I${top_srcdir}/src so += operator to concat CPPFLAGS)
####THISDIR = Pecos/test

AM_CPPFLAGS  = -I${top_srcdir}/src @TEUCHOS_CPPFLAGS@
AM_LDFLAGS   = @TEUCHOS_LDFLAGS@
TEUCHOS_LIBS = -lteuchos

if WITH_GSL
AM_CPPFLAGS += @GSL_CPPFLAGS@
AM_LDFLAGS  += @GSL_LDFLAGS@
GSL_LIBS    =  -lgsl
endif

# Trilinos "model" uses perl scripts to drive tests - for now, simply use 'gmake check'
#### noinst_PROGRAMS = test_pecos_ifft_g test_pecos_ifft_sd

# Tests to run when 'gmake check' is run
TESTS = pecos_ifft_g pecos_ifft_sd
# Tests to build when 'gmake check' is run
check_PROGRAMS = pecos_ifft_g pecos_ifft_sd


# Complete these for each specific test.
pecos_ifft_g_SOURCES = pecos_ifft_g.cpp
pecos_ifft_g_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_g_LDADD = -L${top_builddir} -lpecos \
    -L${top_builddir}/packages/dfftpack -ldfftpack \
    $(TEUCHOS_LIBS) $(GSL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@

# Complete these for each specific test.
pecos_ifft_sd_SOURCES = pecos_ifft_g.cpp
pecos_ifft_sd_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_sd_LDADD = -L${top_builddir} -lpecos \
    -L${top_builddir}/packages/dfftpack -ldfftpack \
    $(TEUCHOS_LIBS) $(GSL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@

# Define this variable if you want to install test header files.
# subpkginclude_HEADERS =
