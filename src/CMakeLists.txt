include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS}
                    ${PROJECT_BINARY_DIR})

# WJB: not sure this is prudent; moving to shared object libs -- link_directories(${Teuchos_LIBRARY_DIRS})

if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
endif(HAVE_SPARSE_GRID)

# The package's pkg_config.h file is located in $PROJECT_BINARY_DIR (see above)
add_definitions("-DHAVE_CONFIG_H")

# WJB: not sure this is prudent; moving to shared object libs -- set(Pecos_libs teuchos lapack)

if(HAVE_FFTW)
  # WJB: not sure this is prudent; moving to shared object libs -- set(Pecos_libs ${Pecos_libs} fftw3)
  include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
endif(HAVE_FFTW)

## Pecos library sources.
file(GLOB Pecos_src *.cpp)

if(HAVE_LHS)
  set(Pecos_src ${Pecos_src} lhs_wrappers.f90)
  set_property(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS ${FIXED_FLAG})
  set(Pecos_libs ${Pecos_libs} lhs)
endif(HAVE_LHS)

add_library(pecos_src SHARED ${Pecos_src})
add_library(pecos_src-static STATIC ${Pecos_src})

## Ensure both the shared and static libraries have the same root name
set_target_properties(pecos_src-static PROPERTIES OUTPUT_NAME "pecos_src")
set_target_properties(pecos_src-static PROPERTIES PREFIX "lib")

target_link_libraries(pecos_src ${Pecos_libs})

