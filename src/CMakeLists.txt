# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

include_directories(${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS}
                    ${Pecos_BINARY_DIR})

#add_definitions("-DHAVE_CONFIG_H")
remove_definitions("-DHAVE_CONFIG_H")

if(HAVE_LHS)
  include_directories(${LHS_BINARY_DIR}/include)
  list(APPEND Pecos_libs lhs mods mod)
endif()

if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
  list(APPEND Pecos_libs sparsegrid)
endif()

if(HAVE_FFT)
  if(HAVE_FFTW)
    list(APPEND Pecos_libs fftw3)
    include_directories(${Pecos_SOURCE_DIR}/packages/fftw/api)
  elseif(HAVE_DFFTPACK)
    list(APPEND Pecos_libs dfftpack)
  endif()
endif()

# Unconditional dependencies.
list(APPEND Pecos_libs ${Teuchos_LIBRARIES}
  ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS}
)

## Pecos library sources.
file(GLOB Pecos_src *.cpp)
file(GLOB Pecos_hdr *.hpp)
list(APPEND Pecos_hdr ${Pecos_BINARY_DIR}/pecos_config.h)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})
install(TARGETS pecos_src EXPORT ${ExportTarget} DESTINATION lib)
install(FILES ${Pecos_hdr} DESTINATION include)
