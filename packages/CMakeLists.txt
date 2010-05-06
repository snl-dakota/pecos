if(NOT USE_SYSTEM_TEUCHOS)
  add_subdirectory(teuchos)
endif()

if(HAVE_FFT)
  add_subdirectory(dfftpack)
#  set(HAVE_FFTW OFF CACHE BOOL "Build the fftw package.  This OPTION has
#  no effect if HAVE_FFT is ON" FORCE)
#  add_subdirectory(fftw)
elseif(HAVE_DFFTPACK)
  add_subdirectory(dfftpack)
elseif(HAVE_FFTW)
  add_subdirectory(fftw)
endif(HAVE_FFT)

if(HAVE_LHS)
  add_subdirectory(LHS)
endif(HAVE_LHS)

if(HAVE_SPARSE_GRID)
  add_subdirectory(VPISparseGrid)
endif(HAVE_SPARSE_GRID)

