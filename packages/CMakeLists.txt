if(NOT USE_SYSTEM_TEUCHOS)
  add_subdirectory(teuchos)
endif()

include(ExternalProject)

if(HAVE_FFT)
  add_subdirectory(dfftpack)
#  set(HAVE_FFTW OFF CACHE BOOL "Build the fftw package.  This OPTION has
#  no effect if HAVE_FFT is ON" FORCE)
#  add_subdirectory(fftw)
elseif(HAVE_DFFTPACK)
  add_subdirectory(dfftpack)
elseif(HAVE_FFTW)
# WJB: use ExternalProj feature (aka heterogeneous bld)  add_subdirectory(fftw)
  set(EXTERNAL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${EXTERNAL_BINARY_DIR}
  )
  find_program(sh sh)
  if(${sh} MATCHES "sh-NOTFOUND")
    message(FATAL_ERROR
      "Cannot build fftw as an external project without sh")
  else()
    #build FFTW as an external project
    externalproject_add(fftw
      DOWNLOAD_COMMAND ""
      CONFIGURE_COMMAND ${sh} ${CMAKE_CURRENT_SOURCE_DIR}/fftw/configure
      BINARY_DIR ${EXTERNAL_BINARY_DIR}
      INSTALL_COMMAND ""
    )
  endif() #sh
endif(HAVE_FFT)

if(HAVE_LHS)
  add_subdirectory(LHS)
endif(HAVE_LHS)

if(HAVE_SPARSE_GRID)
  add_subdirectory(VPISparseGrid)
endif(HAVE_SPARSE_GRID)

