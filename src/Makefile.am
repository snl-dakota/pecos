## Process this file with automake to produce 'Makefile.in'.

AM_CPPFLAGS  = @GSL_CPPFLAGS@ -I${top_builddir}/packages/gsl
AM_CPPFLAGS += @TEUCHOS_CPPFLAGS@ -I${top_builddir}/packages/teuchos/src

## Suppress parallel behavior by gmake; other makes should ignore this.
.NOTPARALLEL:

## Set the Pecos library name.
noinst_LTLIBRARIES = libpecos_src.la

## Pecos library sources.
libpecos_src_la_SOURCES = NatafTransformation.cpp ProbabilityTransformation.cpp

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

