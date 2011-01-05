cmake_minimum_required(VERSION 2.8)

project("Pecos" C CXX Fortran)

include(FortranCInterface)

option(BUILD_SHARED_LIBS "Build Pecos with shared libraries?" OFF)


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


option(BOOST_DISABLE_ASSERTS "Toggle BoostAssert, default is AssertDisabled" ON)
if (BOOST_DISABLE_ASSERTS)
  add_definitions("-DBOOST_DISABLE_ASSERTS")
endif (BOOST_DISABLE_ASSERTS)

# WJB - ToDo: test the (HAVE_FFT == OFF) case
option(HAVE_FFT "Use either the fftw or dfftpack packages" ON)
option(HAVE_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if
  HAVE_FFT is OFF" ON)
option(HAVE_FFTW "Find and use an installed fftw package. This OPTION has no
  effect if HAVE_FFT is OFF" OFF)

option(HAVE_LHS "Build the LHS package." ON)
option(HAVE_SPARSE_GRID "Build the VPISparseGrid package." ON)

#option(USE_SYSTEM_TEUCHOS
#  "If false, Pecos will use the included version of Teuchos."
#  OFF)

#if(USE_SYSTEM_TEUCHOS)
#  # WJB: suspect SYSTEM_TEUCHOS settings specify incorrect paths -- try to get
#  #      a review from someone interested in this case being enabled
#  find_package(Trilinos REQUIRED)
#  ## Kitware version -- include(${Trilinos_DIR}/TrilinosConfig.cmake)
#  include(${Trilinos_DIR}/TeuchosConfig.cmake)
#  ##set(Teuchos_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} CACHE FILEPATH "include dir for Teuchos" FORCE)
#  ##set(Teuchos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS} CACHE FILEPATH "library dir for Teuchos" FORCE)
#else()
#  set(Teuchos_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/packages/teuchos/src
#                           ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/src)
#  set(Teuchos_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/src)
#endif()


#Keep track of what libraries should be included in libpecos
set(SUBLIBS pecos_src)


function(CheckPackage package)
  if(HAVE_${package})
    add_definitions("-DHAVE_${package}")
    string(TOLOWER ${package} package_lower)
    set(SUBLIBS ${SUBLIBS} ${package_lower} PARENT_SCOPE)
  endif()
endfunction(CheckPackage)

#Descend into source subdirectories based on user preference

CheckPackage(LHS)

if(HAVE_FFT)
  CheckPackage(DFFTPACK)

  if(HAVE_FFTW)
    add_definitions("-DHAVE_FFTW")
    include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
    #set(SUBLIBS ${SUBLIBS} fftw3)
    find_package(Fftw3)
  endif(HAVE_FFTW)

endif(HAVE_FFT)

if(HAVE_SPARSE_GRID)
  add_definitions("-DHAVE_SPARSE_GRID")
  set(SUBLIBS ${SUBLIBS} sparsegrid)
endif(HAVE_SPARSE_GRID)


# first probe for system-installed Trilinos
# this will respect Trilinos_DIR if already set
find_package(Trilinos QUIET)
if ( NOT Trilinos_DIR )

  # if no one has configured a local src Teuchos, do so
  # this will respect Teuchos_DIR if already set
  if( NOT Teuchos_DIR )
    set( Teuchos_DIR ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos )
    message( "Setting Teuchos_DIR to ${Teuchos_DIR}" )
    add_subdirectory(packages/teuchos)
  else()
    message( "in ${CMAKE_CURRENT_BINARY_DIR} Teuchos_DIR already set to ${Teuchos_DIR}" )
  endif()
  find_package( Teuchos NO_MODULE REQUIRED )
else()
  message( "Using system trilinos in ${Trilinos_DIR}" )
endif()


add_subdirectory(packages)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# Package-specific configuration options
include_directories(${Pecos_BINARY_DIR})

configure_file("${Pecos_SOURCE_DIR}/cmake/pecos_config.h.in"
  "${Pecos_BINARY_DIR}/pecos_config.h"
  @ONLY
)

add_definitions("-DHAVE_CONFIG_H")


# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
configure_file("${Pecos_SOURCE_DIR}/cmake/empty.cxx.in"
  "${Pecos_BINARY_DIR}/empty.cxx"
  @ONLY
)

add_library(pecos empty.cxx)
target_link_libraries(pecos ${SUBLIBS})

# WJB - ToDo: delete here to bottom once I am done debugging the bldSettings
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt
             "${var} \"${${var}}\"\n")
ENDFOREACH(var ${res})

