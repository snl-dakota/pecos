cmake_minimum_required(VERSION 2.8)

project("Pecos" C CXX Fortran)

include(FortranCInterface)

option(BUILD_SHARED_LIBS "Build Pecos with shared libraries?" ON)
remove_definitions("-DHAVE_CONFIG_H")

# first check for a system blas and lapack
if(NOT DEFINED BLAS_LIBS OR NOT DEFINED LAPACK_LIBS)
  if(NOT DEFINED BLAS_LIBS)
    find_library(BLAS_LIBS blas)
  endif()
  if(NOT DEFINED LAPACK_LIBS)
    find_library(LAPACK_LIBS lapack)
  endif()
  if(NOT BLAS_LIBS OR NOT LAPACK_LIBS)
    # if not a system blas and lapack, then look for a cmake built LAPACK
    # with find_package
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  endif()
endif()

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()


# Mandate a system or user-provided Boost
if(Boost_DIR)
  #message( "in PECOS: Boost_DIR already set to ${Boost_DIR}" )
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
else()
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  endif()
endif(Boost_DIR)


include(CTest)

set( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package_arch
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_tests
  ${CMAKE_MODULE_PATH}
)


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


include(AddExtraTeuchosPaths)

# when building inside Trilinos, the path to Teuchos will already be set
if (NOT BUILD_IN_TRILINOS)
# first probe for system-installed Trilinos
# this will respect Trilinos_DIR if already set
if(NOT Extra_Teuchos_INCLUDE_DIR OR NOT Extra_Teuchos_LIBRARY_DIR)
find_package(Trilinos QUIET)
endif()

if ( NOT Trilinos_DIR )

  # Check a special case in which PATHs to Teuchos headers/lib
  # are specified, independent of a single/common "prefix"
  if(Extra_Teuchos_INCLUDE_DIR AND Extra_Teuchos_LIBRARY_DIR)
    AddExtraTeuchosPaths()
  else()

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
  endif() # Extra_Teuchos_INCLUDE_DIR AND Extra_Teuchos_LIBRARY_DIR
else()
  message( "Using system trilinos in ${Trilinos_DIR}" )
endif() # NOT Trilinos_DIR
endif() # NOT BUILD_IN_TRILINOS
link_directories(${Teuchos_LIBRARY_DIRS})

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

add_subdirectory(packages)

add_subdirectory(src)

option(PECOS_ENABLE_TESTS "Enable Pecos-specific tests?" ON)
if(BUILD_TESTING AND PECOS_ENABLE_TESTS)
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

