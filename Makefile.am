## Process this file with automake to produce 'Makefile.in'.

# AUTOMAKE_OPTIONS = readme-alpha

## Need the config directory for its m4 files.
ACLOCAL_AMFLAGS = -I ./m4

## Source subdirectories to be recursed into during the build.
## Convenience libraries to be built in each of the SUBDIRS.
## Top-level, libpecos.la must be built PRIOR to the linkage of test executables
SUBDIRS = packages src . test
SUBLIBS = src/libpecos_src.la

if WITH_FFT
SUBLIBS += packages/dfftpack/libdfftpack.la packages/fftw/libfftw.la
else
if WITH_DFFTPACK
SUBLIBS += packages/dfftpack/libdfftpack.la
endif
if WITH_FFTW
SUBLIBS += packages/fftw/libfftw3.la
endif
endif
if WITH_LHS
SUBLIBS += packages/LHS/liblhs.la
endif
if WITH_VPI_SPARSE_GRID
SUBLIBS += packages/VPISparseGrid/src/libvpisparsegrid.la
endif

## Set the Pecos library name.
lib_LTLIBRARIES = libpecos.la

## Pecos headers that need to be distributed
pkginclude_HEADERS = pecos_config.h				\
		     src/pecos_data_types.hpp			\
		     src/pecos_global_defs.hpp			\
		     src/pecos_stat_util.hpp			\
		     src/DataTransformation.hpp			\
		     src/InverseTransformation.hpp		\
		     src/FourierInverseTransformation.hpp	\
		     src/KarhunenLoeveInverseTransformation.hpp	\
		     src/SamplingInverseTransformation.hpp	\
		     src/ProbabilityTransformation.hpp		\
		     src/NatafTransformation.hpp		\
		     src/LHSDriver.hpp				\
		     src/BoostRNG_Monostate.hpp
#		     src/BasisFunction.hpp

## No sources in this directory!
libpecos_la_SOURCES = 

## Combine the convenience libraries into the final product.
libpecos_la_LIBADD = $(SUBLIBS)

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = aclocal.m4 config.guess config.sub configure \
		       pecos_config.h.in install-sh ltmain.sh missing Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
#EXTRA_DIST = configure.ac README.pecos

