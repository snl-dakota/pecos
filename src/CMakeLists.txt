# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

add_definitions("-DHAVE_CONFIG_H")

include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS}
                    ${PROJECT_BINARY_DIR})

if(HAVE_LHS)
  list(APPEND Pecos_libs lhs)
endif()

if(HAVE_SPARSE_GRID)
  list(APPEND Pecos_libs sparsegrid)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
endif()

if(HAVE_FFT)
  if(HAVE_FFTW)
    list(APPEND Pecos_libs fftw3)
    include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/api)
  elseif(HAVE_DFFTPACK)
    list(APPEND Pecos_libs dfftpack)
  endif()
endif()

# Unconditional dependencies.
list(APPEND Pecos_libs
  ${Teuchos_LIBRARIES}
  ${LAPACK_LIBS}
  ${BLAS_LIBS}
  ${FLIBS}
  ${FCLIBS}
  )

## Pecos library sources.
file(GLOB Pecos_src *.cpp)
file(GLOB Pecos_hdr *.hpp)
list(APPEND Pecos_hdr ${Pecos_BINARY_DIR}/pecos_config.h)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})
install(TARGETS pecos_src DESTINATION lib)
install(FILES ${Pecos_hdr} DESTINATION include)
