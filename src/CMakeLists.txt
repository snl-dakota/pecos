# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

add_definitions("-DHAVE_CONFIG_H")

include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS}
                    ${PROJECT_BINARY_DIR})


if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
endif(HAVE_SPARSE_GRID)

if(HAVE_FFTW)
  #set(Pecos_libs ${Pecos_libs} fftw3)
  include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
endif(HAVE_FFTW)

## Pecos library sources.
file(GLOB Pecos_src *.cpp)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})
