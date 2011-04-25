cmake_minimum_required(VERSION 2.8)

project("Pecos" C CXX Fortran)

include(FortranCInterface)

option(BUILD_SHARED_LIBS "Build Pecos with shared libraries?" ON)
remove_definitions("-DHAVE_CONFIG_H")


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

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()


# Mandate a system or user-provided Boost
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
if(HAVE_LHS)
  list(APPEND SUBLIBS mods mod)
endif(HAVE_LHS)

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


# when building inside Trilinos, the path to Teuchos will already be set
if (NOT BUILD_IN_TRILINOS)
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
    set(Teuchos_LIBRARY_DIRS "")
  else()
    message( "in ${CMAKE_CURRENT_BINARY_DIR} Teuchos_DIR already set to ${Teuchos_DIR}" )
  endif()
  set(Teuchos_LIBRARIES teuchos)
  find_package( Teuchos NO_MODULE REQUIRED )
else()
  message( "Using system trilinos in ${Trilinos_DIR}" )
endif()
endif()
link_directories(${Teuchos_LIBRARY_DIRS})

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

add_subdirectory(packages)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()


# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
configure_file("${Pecos_SOURCE_DIR}/cmake/empty.cxx.in"
  "${Pecos_BINARY_DIR}/empty.cxx"
  @ONLY
)

add_library(pecos empty.cxx)
target_link_libraries(pecos ${SUBLIBS})
install(TARGETS pecos EXPORT ${ExportTarget} DESTINATION lib)

# WJB - ToDo: delete here to bottom once I am done debugging the bldSettings
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt
             "${var} \"${${var}}\"\n")
ENDFOREACH(var ${res})

set(Pecos_LINK_LIBRARIES ${SUBLIBS} CACHE INTERNAL "Pecos link libraries")

