# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

#the following configuration file is necessary to build shared libs on windows
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PecosDLLExportMacro.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/PecosDLLExportMacro.h)
  
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#keep track of the libraries we need to link against
set(Pecos_libs "teuchos")

add_definitions("-DHAVE_CONFIG_H")

include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS}
                    ${PROJECT_BINARY_DIR}
					${CMAKE_CURRENT_BINARY_DIR})


if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
  include_directories(${PROJECT_BINARY_DIR}/packages/VPISparseGrid/src)
  list(APPEND Pecos_libs "sparsegrid")
endif()

if(HAVE_FFTW)
  #set(Pecos_libs ${Pecos_libs} fftw3)
  include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
endif()

if(HAVE_LHS)
  list(APPEND Pecos_libs "lhs")
endif()

if(HAVE_DFFTPACK)
  list(APPEND Pecos_libs "dfftpack")
endif()

## Pecos library sources.
file(GLOB Pecos_src *.cpp)

add_library(pecos ${Pecos_src})
target_link_libraries(pecos ${Pecos_libs})
