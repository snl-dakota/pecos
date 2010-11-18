## Process this file with automake to produce 'Makefile.in'.

## Initialize the SUBDIRS variable; it'll be added to below.
SUBDIRS = 

if BUILD_TEUCHOS
SUBDIRS += teuchos
endif

if BUILD_BOOST_TPL
SUBDIRS += boost
endif

## Package subdirectories to recurse into.
if WITH_FFT
SUBDIRS += dfftpack fftw
else
if WITH_DFFTPACK
SUBDIRS += dfftpack
endif
if WITH_FFTW
SUBDIRS += fftw
endif
endif

if WITH_LHS
SUBDIRS += LHS
endif

if WITH_SPARSE_GRID
SUBDIRS += VPISparseGrid
endif

## Conditionally configured subdirectories should _not_ be recursed into
## during 'make distclean' or 'make maintainer-clean'.
DIST_SUBDIRS = $(SUBDIRS)

## These files should be deleted on 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

