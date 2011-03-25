## Process this file with automake to produce 'Makefile.in'.

if WITH_DEBUG_ENABLED
AM_CXXFLAGS = @DEBUGFLAGS@
endif

AM_CPPFLAGS = @TEUCHOS_CPPFLAGS@ @BOOST_CPPFLAGS@

if WITH_FFTW
AM_CPPFLAGS += @FFTW_CPPFLAGS@ -I${top_srcdir}/packages/fftw/api
endif

if WITH_SPARSE_GRID
AM_CPPFLAGS += -I${top_srcdir}/packages/VPISparseGrid/src
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
			  BasisApproximation.cpp			\
			  PolynomialApproximation.cpp			\
			  InterpPolyApproximation.cpp			\
			  LagrangeInterpPolyApproximation.cpp		\
			  HermiteInterpPolyApproximation.cpp		\
			  OrthogPolyApproximation.cpp			\
			  BasisPolynomial.cpp				\
			  InterpolationPolynomial.cpp			\
			  LagrangeInterpPolynomial.cpp			\
			  PiecewiseInterpPolynomial.cpp			\
			  HierarchicalLinearBasis.cpp			\
			  OrthogonalPolynomial.cpp			\
			  ChebyshevOrthogPolynomial.cpp			\
			  GenLaguerreOrthogPolynomial.cpp		\
			  HermiteOrthogPolynomial.cpp			\
			  JacobiOrthogPolynomial.cpp			\
			  LaguerreOrthogPolynomial.cpp			\
			  LegendreOrthogPolynomial.cpp			\
			  NumericGenOrthogPolynomial.cpp		\
                          LHSDriver.cpp					\
			  IntegrationDriver.cpp				\
			  CubatureDriver.cpp				\
			  SparseGridDriver.cpp				\
			  TensorProductDriver.cpp			\
			  pecos_stat_util.cpp

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
