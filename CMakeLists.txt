cmake_minimum_required(VERSION 2.8)

PROJECT("Pecos" C CXX Fortran)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF() 

INCLUDE(CTest)

SET( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package_arch
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_tests
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/packages/boost)

CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)

# Package-specific configuration options
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/pecos_config.h.in"
  "${PROJECT_BINARY_DIR}/pecos_config.h"
)

LIST(APPEND CMAKE_CXX_FLAGS "-DHAVE_CONFIG_H -DHAVE_SPARSE_GRID")

SET(Pecos_ENABLE_C ON)
SET(Pecos_ENABLE_CXX ON)
SET(Pecos_ENABLE_Fortran ON)

#INCLUDE(TrilinosFortranMangling)

SET(Pecos_ENABLE_DEBUG OFF)
SET(Pecos_ENABLE_EXPLICIT_INSTANTIATION OFF)
SET(Pecos_ENABLE_PERFORMANCE_TESTS OFF)

SUBDIRS(packages src test)

# WJB - ToDo: upgrade to latest FFTW prior to CMakeifying
#OPTION(WITH_FFT "Build the fftw and dfftpack packages" ON)
OPTION(WITH_FFT "Build the dfftpack packages" ON)
OPTION(WITH_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if WITH_FFT is ON" ON)
OPTION(WITH_FFTW "Build the fftw package.  This OPTION has no effect if WITH_FFT is ON" OFF)
OPTION(WITH_LHS "Build the LHS package." ON)

#Keep track of what libraries should be included in libpecos
SET(SUBLIBS pecos_src)

#Descend into source subdirectories based on user preference
IF(WITH_FFT)
#  SET(SUBLIBS ${SUBLIBS} dfftpack fftw3)
  SET(SUBLIBS ${SUBLIBS} dfftpack)
  ADD_DEFINITIONS("-DHAVE_DFFTPACK")
ELSEIF(WITH_DFFTPACK)
  SET(SUBLIBS ${SUBLIBS} dfftpack)
  ADD_DEFINITIONS("-DHAVE_DFFTPACK")
ELSEIF(WITH_FFTW)
  SET(SUBLIBS ${SUBLIBS} fftw3)
ENDIF(WITH_FFT)

# WJB - ToDo: get a consensus on whether 'HAVE_LHS' should be defined in
# pecos_config.h or on the compile line
IF(WITH_LHS)
  SET(SUBLIBS ${SUBLIBS} lhs)
  ADD_DEFINITIONS("-DHAVE_LHS")
ENDIF(WITH_LHS)

# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
ADD_LIBRARY(pecos empty.cxx)
TARGET_LINK_LIBRARIES(pecos ${SUBLIBS})

#not sure what's up with this, install rules or something???
## Pecos headers that need to be distributed
#pkginclude_HEADERS = pecos_config.h				\
#		     src/pecos_data_types.hpp			\
#		     src/pecos_global_defs.hpp			\
#		     src/pecos_stat_util.hpp			\
#		     src/DataTransformation.hpp			\
#		     src/InverseTransformation.hpp		\
#		     src/FourierInverseTransformation.hpp	\
#		     src/KarhunenLoeveInverseTransformation.hpp	\
#		     src/SamplingInverseTransformation.hpp	\
#		     src/ProbabilityTransformation.hpp		\
#		     src/NatafTransformation.hpp		\
#		     src/LHSDriver.hpp				\
#		     src/BoostRNG_Monostate.hpp
##		     src/BasisFunction.hpp

