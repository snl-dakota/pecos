## Process this file with automake to produce 'Makefile.in'.

# WJB - ToDo ASAP: proper treatment of path to "external" GSL (possibly in /usr/include)
AM_CPPFLAGS  = -I${top_builddir}/packages/gsl @GSL_CPPFLAGS@

if WITH_ALT_EXTERNAL_TEUCHOS
AM_CPPFLAGS += @TEUCHOS_CPPFLAGS@
else
AM_CPPFLAGS += -I${top_builddir}/packages/teuchos/src
endif

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
			  LHSDriver.cpp
#			  BasisFunction.cpp

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
