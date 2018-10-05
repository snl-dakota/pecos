#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

cmake_minimum_required(VERSION 3.1)

set( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

# Pecos requires C++11
include(DakotaCxxOptions)
dakota_cxx_standard()

project("Pecos" C CXX Fortran)

# --- CMake modules ---
include(FortranCInterface)
include(CTest)


# --- Global options ---

option(BUILD_SHARED_LIBS "Build Pecos with shared libraries?" ON)

remove_definitions("-DHAVE_CONFIG_H")

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")


# --- External TPLs ---

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

# Mandate a system or user-provided Boost
if(Boost_DIR)
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
else()
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  endif()
endif(Boost_DIR)

if (NOT BUILD_IN_TRILINOS)

  # First probe for system-installed Trilinos, respecting Trilinos_DIR if set
  find_package(Trilinos QUIET)

  if(Trilinos_FOUND)
    if(Trilinos_DIR)
      message(STATUS
	"Pecos using previously specified Trilinos in ${Trilinos_DIR}")
    else()
      message(STATUS "Pecos using external Trilinos")
    endif()
  else()

    # If no parent project configured Teuchos, do so, using Teuchos_DIR if set
    if(Teuchos_DIR)
      message(STATUS
	"Pecos using previously specified Teuchos in ${Teuchos_DIR}")
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/packages/trilinos)

      set(Teuchos_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/packages/trilinos/packages/teuchos)
      set(Trilinos_ENABLE_Teuchos ON CACHE BOOL
        "Pecos enabling Trilinos Teuchos" FORCE)

      # Map key Pecos variables to TriBITS variables
      set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
      set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )
      # Pecos doesn't use any Teuchos MPI features; may want to force off
      #set( TPL_ENABLE_MPI ${DAKOTA_HAVE_MPI} )
      # This doesn't do as name implies; setting OFF doesn't generate Config.cmake
      # at all; doesn't just control whether installed!  Want Config.cmake in build
      #        set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL
      #	  "Pecos is the top-level package; don't write Trilinos config files")

      message(STATUS "Pecos setting Teuchos_DIR to ${Teuchos_DIR}")

      # This mirrors the Trilinos release process setting, so Pecos
      # can work with Trilinos master, when needed.  Specifically it
      # will suppress errors about missing packages.
      set(Trilinos_ENABLE_DEVELOPMENT_MODE OFF CACHE BOOL
	"Pecos disabling Trilinos development mode")
      # This would be lighter weight, but might result in strong warnings
      ##set(Trilinos_ASSERT_MISSING_PACKAGES OFF CACHE BOOL "Dakota being lazy")

      set(Trilinos_GENERATE_REPO_VERSION_FILE OFF CACHE BOOL
        "Pecos disabling generation of TrilinosRepoVersion.txt")

      # TODO: Propagate to top-level Dakota configure
      set(Teuchos_ENABLE_COMPLEX TRUE CACHE BOOL
	"Pecos requires some complex blas wrappers")

      add_subdirectory(packages/trilinos)

    else()
      message(WARNING
	"Pecos requires teuchos, but could not find local nor parent copy!")
    endif() # Teuchos_DIR

    # Additional setting to prevent multiple targets with the same name
    set(Trilinos_TARGETS_IMPORTED 1)

    find_package( Teuchos NO_MODULE REQUIRED )

  endif() # Trilinos_DIR

endif() # NOT BUILD_IN_TRILINOS

link_directories(${Teuchos_LIBRARY_DIRS})



# --- Options for Pecos components
##option(PECOS_ENABLE_SAMPLING "Enable sampling module" ON)
option(PECOS_ENABLE_ORTHOGPOLY "Enable core pecos module" ON)
option(PECOS_ENABLE_SURROGATES "Enable surrogates module" ON)
option(PECOS_ENABLE_PYTHONAPI "Enable Python wrappers" OFF)
option(PECOS_ENABLE_TESTS "Enable Pecos-specific tests?" ON)

# --- Potential locally included packages ---

option(HAVE_FFT "Use either the fftw or dfftpack packages" ON)
option(HAVE_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if
  HAVE_FFT is OFF" ON)
option(HAVE_FFTW "Find and use an installed fftw package. This OPTION has no
  effect if HAVE_FFT is OFF" OFF)

option(HAVE_LHS "Build the LHS package." ON)
option(HAVE_SPARSE_GRID "Build the VPISparseGrid package." ON)

add_subdirectory(packages)


# --- Dependency configuration based on previously defined packages ---

# Keep track of what libraries should be included in libpecos
set(SUBLIBS pecos_util)

if(PECOS_ENABLE_ORTHOGPOLY)
  list(APPEND SUBLIBS pecos_src)
endif()

function(CheckPackage package)
  if(HAVE_${package})
    add_definitions("-DHAVE_${package}")
    string(TOLOWER ${package} package_lower)
    set(SUBLIBS ${SUBLIBS} ${package_lower} PARENT_SCOPE)
  endif()
endfunction(CheckPackage)

CheckPackage(LHS)
if(HAVE_LHS)
  list(APPEND SUBLIBS lhs_mods lhs_mod)
endif(HAVE_LHS)

if(HAVE_FFT)
  CheckPackage(DFFTPACK)

  if(HAVE_FFTW)
    # BMA: I couldn't get this to work with the ExternalProject; commenting:
    #find_package(Fftw3)
    list(APPEND SUBLIBS fftw3)
    link_directories("${FFTW_BINARY_DIR}/.libs")
  endif(HAVE_FFTW)

endif(HAVE_FFT)

if(HAVE_SPARSE_GRID)
  add_definitions("-DHAVE_SPARSE_GRID")
  set(SUBLIBS ${SUBLIBS} sparsegrid)
endif(HAVE_SPARSE_GRID)



# Foundational utilities
add_subdirectory(util)

# Core Pecos
if(PECOS_ENABLE_ORTHOGPOLY)
  add_subdirectory(src)
endif()

if(BUILD_TESTING AND PECOS_ENABLE_TESTS AND PECOS_ENABLE_ORTHOGPOLY)
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

set(Pecos_LINK_LIBRARIES ${SUBLIBS} CACHE INTERNAL "Pecos link libraries")

# Modular surrogate model library
# TODO: optionally depences on Pecos core (need tests with and without)
if(PECOS_ENABLE_SURROGATES)
  add_subdirectory(surrogates)
endif()

# Python wrappers are optional; this approach is having CMake invoke
# distutils setup.py. Must add this subdirectory after all other core
# libraries. TODO: check for required dependencies.
if(PECOS_ENABLE_PYTHONAPI)
  find_package(PythonInterp REQUIRED)
  add_subdirectory(python)
endif()
