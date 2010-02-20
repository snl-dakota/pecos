cmake_minimum_required(VERSION 2.6)

PROJECT("Pecos" C CXX Fortran)

IF(WIN32)
  FIND_PACKAGE(CLAPACK REQUIRED)
ELSE()
  FIND_PACKAGE(BLAS REQUIRED)
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

SET(Pecos_ENABLE_C ON)
SET(Pecos_ENABLE_CXX ON)
SET(Pecos_ENABLE_Fortran ON)

SUBDIRS(packages src test)

#INCLUDE(TrilinosFortranMangling)

# Package-specific configuration options
SET(Pecos_ENABLE_DEBUG OFF)
SET(Pecos_ENABLE_EXPLICIT_INSTANTIATION OFF)
SET(Pecos_ENABLE_PERFORMANCE_TESTS OFF)

SET(BOOST_DISABLE_ASSERTS 1)

# WJB - ToDo: upgrade to latest FFTW prior to CMakeifying
#OPTION(WITH_FFT "Build the fftw and dfftpack packages" ON)
OPTION(WITH_FFT "Build the dfftpack packages" ON)
OPTION(HAVE_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if WITH_FFT is ON" ON)
OPTION(WITH_FFTW "Build the fftw package.  This OPTION has no effect if WITH_FFT is ON" OFF)

# WJB: why 2 vars when 1 will do?? OPTION(WITH_LHS "Build the LHS package." ON)
OPTION(HAVE_LHS "Build the LHS package." ON)
OPTION(HAVE_SPARSE_GRID "Build the VPISparseGrid package." ON)

#Keep track of what libraries should be included in libpecos
SET(SUBLIBS pecos_src)

#Descend into source subdirectories based on user preference
IF(WITH_FFT)
#  SET(SUBLIBS ${SUBLIBS} dfftpack fftw3)
  SET(SUBLIBS ${SUBLIBS} dfftpack)
ELSEIF(HAVE_DFFTPACK)
  SET(SUBLIBS ${SUBLIBS} dfftpack)
ELSEIF(WITH_FFTW)
  SET(SUBLIBS ${SUBLIBS} fftw3)
ENDIF(WITH_FFT)

IF(HAVE_LHS)
  SET(SUBLIBS ${SUBLIBS} lhs)
ENDIF(HAVE_LHS)

IF(HAVE_SPARSE_GRID)
  SET(SUBLIBS ${SUBLIBS} sparsegrid)
ENDIF(HAVE_SPARSE_GRID)


# Package-specific configuration options
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/cmake/pecos_config.h.in"
  "${PROJECT_BINARY_DIR}/pecos_config.h"
  @ONLY
)

ADD_DEFINITIONS("-DHAVE_CONFIG_H")


# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
ADD_LIBRARY(pecos empty.cxx)
TARGET_LINK_LIBRARIES(pecos ${SUBLIBS})

# WJB - ToDo: delete here to bottom once I am done debugging the build settings
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt
             "${var} \"${${var}}\"\n")
ENDFOREACH(var ${res})

