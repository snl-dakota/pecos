## Process this file with automake to produce 'Makefile.in'.

## Initialize the SUBDIRS variable; it'll be added to below.
if WITH_ALT_EXTERNAL_TEUCHOS
SUBDIRS = 
else
SUBDIRS = teuchos
endif

## Package subdirectories to recurse into.
if WITH_GSL
SUBDIRS += gsl
endif

#if WITH_FFT
#SUBDIRS += dfftpack fftw
#else
#  if WITH_DFFTPACK
#  SUBDIRS += dfftpack
#  endif
#  if WITH_FFTW
#  SUBDIRS += fftw
#  endif
#endif

## Conditionally configured subdirectories should _not_ be recursed into
## during 'make distclean' or 'make maintainer-clean'.
DIST_SUBDIRS = $(SUBDIRS)

## These files should be deleted on 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

