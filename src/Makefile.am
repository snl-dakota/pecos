## Process this file with automake to produce 'Makefile.in'.

AM_CPPFLAGS = @TEUCHOS_CPPFLAGS@

if WITH_BOOST
AM_CPPFLAGS += @BOOST_CPPFLAGS@ -I$(top_srcdir)/packages
endif

if WITH_GSL
AM_CPPFLAGS += @GSL_CPPFLAGS@
endif

if WITH_FFTW
AM_CPPFLAGS += @FFTW_CPPFLAGS@ -I$(top_srcdir)/packages/fftw/api
endif

if WITH_LHS
LHS_S         = lhs_wrappers.f90
endif

## HACK NOTICE: Libtool 1.5.X doesn't provide tagging support for F90 sources.
## Until it does (Libtool 2.0.X), we need to use the F77 tag rules.
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) \
	      $(FCFLAGS) $(FCFLAGS_f90) $(MODULES)

## Set the Pecos library name.
noinst_LTLIBRARIES = libpecos_src.la

## Pecos library sources.
libpecos_src_la_SOURCES = DataTransformation.cpp			\
			  InverseTransformation.cpp			\
			  FourierInverseTransformation.cpp		\
			  KarhunenLoeveInverseTransformation.cpp	\
			  SamplingInverseTransformation.cpp		\
			  ProbabilityTransformation.cpp			\
			  NatafTransformation.cpp			\
                          Random.cpp                                    \
			  LHSDriver.cpp $(LHS_S)
#			  BasisFunction.cpp

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
