#include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS})
link_directories(${Teuchos_LIBRARY_DIRS})

if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
endif(HAVE_SPARSE_GRID)

set(Pecos_libs teuchos lapack)

if(WITH_FFTW)
  set(Pecos_libs ${Pecos_libs} fftw3)
  include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
endif(WITH_FFTW)

## Pecos library sources.
file(GLOB Pecos_src *.cpp)

if(WITH_LHS)
  set(Pecos_src ${Pecos_src} lhs_wrappers.f90)
  set_property(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS ${FIXED_FLAG})
  set(Pecos_libs ${Pecos_libs} lhs)
endif(WITH_LHS)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})

