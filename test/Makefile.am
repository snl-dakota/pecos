# WJB - ToDo: look into the possibility of "sharing" across all Pecos Makefiles
####include $(top_srcdir)/Pecos.am.mk

# This definition is needed to make full use of Pecos.am.mk, e.g init AM_CPPFLAGS (set to -I../src so += operator works)
####THISDIR = Pecos/test

AM_CPPFLAGS = -I../src @TEUCHOS_CPPFLAGS@

# Tests to run when 'gmake check' is run
TESTS = pecos_ifft_g pecos_ifft_sd
# Tests to build when 'gmake check' is run
check_PROGRAMS = pecos_ifft_g pecos_ifft_sd

# Complete these for each specific test.
pecos_ifft_g_SOURCES = pecos_ifft_g.cpp
#pecos_ifft_g_DEPENDENCIES = ../packages/libgsl.la
#pecos_ifft_g_LDADD = -L.. -lgsl
pecos_ifft_g_LDADD = -L.. -lpecos -L../packages/dfftpack -ldfftpack \
    @TEUCHOS_LDFLAGS@ -lteuchos

# Complete these for each specific test.
pecos_ifft_sd_SOURCES = pecos_ifft_g.cpp
#pecos_ifft_sd_DEPENDENCIES = ../packages/libgsl.la
#pecos_ifft_sd_LDADD = -L.. -lgsl
pecos_ifft_sd_LDADD = -L.. -lpecos -L../packages/dfftpack -ldfftpack \
    @TEUCHOS_LDFLAGS@ -lteuchos

# Define this variable if you want to install test header files.
# subpkginclude_HEADERS =
