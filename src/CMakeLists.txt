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
file(GLOB Pecos_hdr *.hpp)
list(APPEND Pecos_hdr ${Pecos_BINARY_DIR}/pecos_config.h)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})
install(TARGETS pecos_src DESTINATION lib)
install(FILES ${Pecos_hdr} DESTINATION include)
