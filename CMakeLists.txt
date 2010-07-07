cmake_minimum_required(VERSION 2.8)

project("Pecos" C CXX Fortran)

# WJB: dependency on Trilinos macros seems undesirable -- #include(TrilinosFortranMangling)
include(FortranCInterface)

option(BUILD_STATIC "build static libraries" OFF)


# first check for a system blas and lapack
find_library(BLAS_LIBS blas)
find_library(LAPACK_LIBS lapack)
if(NOT BLAS_LIBS OR NOT LAPACK_LIBS)
  # if not a system blas and lapack, then
  # look for a cmake built LAPACK with
  # find_package
  find_package(LAPACK REQUIRED NO_MODULE)
  set(BLAS_LIBS blas)
  set(LAPACK_LIBS lapack)
endif()

# hack to be fixed when CMake does it automatically
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(FIXED_FLAG "-fixed")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(FIXED_FLAG "-ffixed-form")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set(FIXED_FLAG "-Mfixed")
endif()

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

find_package(Boost 1.40 REQUIRED)

include(CTest)

set( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package_arch
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_tests
  ${CMAKE_MODULE_PATH}
)


# should not be hard coded
set(Pecos_ENABLE_DEBUG ON)
set(Pecos_ENABLE_EXPLICIT_INSTANTIATION OFF)
set(Pecos_ENABLE_PERFORMANCE_TESTS OFF)

set(BOOST_DISABLE_ASSERTS 1)

# WJB: see conditional further down -- find_package(Fftw3)

# WJB - ToDo: upgrade to latest FFTW prior to FULL CMakeification
#option(HAVE_FFT "Build the fftw and dfftpack packages" ON)
option(HAVE_FFT "FIND an installed fftw pkg and BUILD the dfftpack pkg in the
distro" ON)
option(HAVE_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if
  WITH_FFT is ON" ON)
# WJB - ToDo: ensure that the default setting for HAVE_FFTW is OFF, even when
# HAVE_FFT is ON
option(HAVE_FFTW "FIND and installed fftw package. This OPTION has no effect
if HAVE_FFT is ON" OFF)

option(HAVE_LHS "Build the LHS package." ON)
option(HAVE_SPARSE_GRID "Build the VPISparseGrid package." ON)

option(USE_SYSTEM_TEUCHOS
  "If false, Pecos will use the included version of Teuchos."
  OFF)

if(USE_SYSTEM_TEUCHOS)
  # WJB: suspect SYSTEM_TEUCHOS settings specify incorrect paths -- try to get
  #      a review from someone interested in this case being enabled
  find_package(Trilinos REQUIRED)
  ## Kitware version -- include(${Trilinos_DIR}/TrilinosConfig.cmake)
  include(${Trilinos_DIR}/TeuchosConfig.cmake)
  ##set(Teuchos_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} CACHE FILEPATH "include dir for Teuchos" FORCE)
  ##set(Teuchos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS} CACHE FILEPATH "library dir for Teuchos" FORCE)
else()
  set(Teuchos_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/packages/teuchos/src
                           ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/src)
  set(Teuchos_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/src)
endif()


#Keep track of what libraries should be included in libpecos
set(SUBLIBS pecos_src)

#Descend into source subdirectories based on user preference
if(HAVE_FFT)
  #set(SUBLIBS ${SUBLIBS} dfftpack fftw3)
  set(SUBLIBS ${SUBLIBS} dfftpack)
elseif(HAVE_DFFTPACK)
  set(SUBLIBS ${SUBLIBS} dfftpack)
elseif(HAVE_FFTW)
  #set(SUBLIBS ${SUBLIBS} fftw3)
  find_package(Fftw3)
endif(HAVE_FFT)

if(HAVE_LHS)
  set(SUBLIBS ${SUBLIBS} lhs)
endif(HAVE_LHS)

if(HAVE_SPARSE_GRID)
  set(SUBLIBS ${SUBLIBS} sparsegrid)
endif(HAVE_SPARSE_GRID)

add_subdirectory(packages)
add_subdirectory(src)
add_subdirectory(test)

# Package-specific configuration options
#include_directories(${Pecos_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR})

configure_file("${PROJECT_SOURCE_DIR}/cmake/pecos_config.h.in"
  "${PROJECT_BINARY_DIR}/pecos_config.h"
  @ONLY
)

add_definitions("-DHAVE_CONFIG_H")


# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
add_library(pecos SHARED empty.cxx)
add_library(pecos-static STATIC empty.cxx)

## Ensure both the shared and static libraries have the same root name
set_target_properties(pecos-static PROPERTIES OUTPUT_NAME "pecos")
set_target_properties(pecos-static PROPERTIES PREFIX "lib")

target_link_libraries(pecos ${SUBLIBS})


# WJB - ToDo: delete here to bottom once I am done debugging the build
# settings
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt
             "${var} \"${${var}}\"\n")
ENDFOREACH(var ${res})

