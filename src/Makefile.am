## Process this file with automake to produce 'Makefile.in'.

## Need the m4 directory for non-native macro tests. (WJB: even in this subdir?)
#ACLOCAL_AMFLAGS = -I ../../../m4

AM_CPPFLAGS  = @GSL_CPPFLAGS@ -I${top_builddir}/packages/gsl
AM_CPPFLAGS += @TEUCHOS_CPPFLAGS@ -I${top_builddir}/packages/teuchos/src

## Suppress parallel behavior by gmake; other makes should ignore this.
.NOTPARALLEL:

## Set the Pecos (WJB: note in a subdir!)library name.
noinst_LTLIBRARIES = libpecos_src.la

## Pecos library sources.
libpecos_src_la_SOURCES = NatafTransformation.cpp Transformation.cpp

## WJB: remove specification of headers if better handled at the top of Pecos
#include_HEADERS = Transformation.hpp pecos_data_types.hpp

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

