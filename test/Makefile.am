## Process this file with automake to produce 'Makefile.in'.

if WITH_DEBUG_ENABLED
AM_CXXFLAGS = @DEBUGFLAGS@
endif

AM_CPPFLAGS  = -I${top_srcdir}/src @TEUCHOS_CPPFLAGS@ @BOOST_CPPFLAGS@
TPL_LDFLAGS  = @TEUCHOS_LDFLAGS@
TPL_LIBS     = -lteuchos

if WITH_LHS
TPL_LDFLAGS += @LHS_LDFLAGS@
TPL_LIBS    += -llhs
endif

if WITH_DFFTPACK
TPL_LDFLAGS += @DFFTPACK_LDFLAGS@
TPL_LIBS    += -ldfftpack
endif

if WITH_FFTW
TPL_LDFLAGS += @FFTW_LDFLAGS@
TPL_LIBS    += -lfftw3
endif

if WITH_SPARSE_GRID
AM_CPPFLAGS += -I${top_srcdir}/packages/VPISparseGrid/src
TPL_LDFLAGS += -L${top_builddir}/packages/VPISparseGrid/src
TPL_LIBS    += -lsparsegrid
endif


# Tests to run when 'gmake check' is run

TESTS = pecos_int_driver pecos_lhs_driver pecos_ifft_g pecos_ifft_sd \
	boost_test_dist boost_test_rng

# Programs to link when 'gmake' is run
noinst_PROGRAMS = pecos_int_driver pecos_lhs_driver pecos_ifft_g pecos_ifft_sd \
		  boost_test_dist boost_test_rng

#if HAVE_CPPUNIT
#TESTS += pecos_hierarchical_driver
#noinst_PROGRAMS += pecos_hierarchical_driver
#endif

# Complete these for each specific test.
pecos_int_driver_SOURCES = pecos_int_driver.cpp
pecos_int_driver_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_int_driver_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

#pecos_hierarchical_test_SOURCES = pecos_hierarchical_driver.cpp HierarchicalLinearBasisTestClass.cpp
#pecos_hierarchical_test_DEPENDENCIES = ${top_builddir}/libpecos.la
#pecos_hierarchical_test_LDADD = -L${top_builddir} -lpecos -L${CPPUNIT_LIBS} -lstdc++ -lcppunit -ldl \
#    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

pecos_lhs_driver_SOURCES = pecos_lhs_driver.cpp
pecos_lhs_driver_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_lhs_driver_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

pecos_ifft_g_SOURCES = pecos_ifft_g.cpp
pecos_ifft_g_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_g_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

pecos_ifft_sd_SOURCES = pecos_ifft_sd.cpp
pecos_ifft_sd_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_ifft_sd_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

boost_test_dist_SOURCES = boost_test_dist.cpp
boost_test_dist_DEPENDENCIES = ${top_builddir}/libpecos.la
boost_test_dist_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@

boost_test_rng_SOURCES = boost_test_rng.cpp
boost_test_rng_DEPENDENCIES = ${top_builddir}/libpecos.la
boost_test_rng_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@


check_PROGRAMS = $(TESTS)
pecos_hierarchical_driver_SOURCES = pecos_hierarchical_driver.cpp HierarchicalLinearBasisTestClass.cpp
pecos_hierarchical_driver_DEPENDENCIES = ${top_builddir}/libpecos.la
pecos_hierarchical_driver_CXXFLAGS = $(CPPUNIT_CFLAGS)
pecos_hierarchical_driver_LDFLAGS = $(CPPUNIT_LIBS)
pecos_hierarchical_driver_LDADD = -L${top_builddir} -lpecos \
    $(TPL_LDFLAGS) $(TPL_LIBS) @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @FCLIBS@


# Define this variable if you want to install test header files.
# subpkginclude_HEADERS =
