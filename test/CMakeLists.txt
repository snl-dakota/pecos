set(TPL_LIBS teuchos)

include_directories(${PROJECT_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS})

link_directories(${Teuchos_LIBRARY_DIRS})

if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
  set(TPL_LIBS ${TPL_LIBS} sparsegrid)
endif(HAVE_SPARSE_GRID)

if(HAVE_LHS)
  set(TPL_LIBS ${TPL_LIBS} lhs)
endif(HAVE_LHS)

if(HAVE_DFFTPACK)
  set(TPL_LIBS ${TPL_LIBS} dfftpack)
endif(HAVE_DFFTPACK)

if(HAVE_FFTW)
  set(TPL_LIBS ${TPL_LIBS} fftw3)
endif(HAVE_FFTW)

# Complete these for each specific test.
add_executable(pecos_ifft_g pecos_ifft_g.cpp)
target_link_libraries(pecos_ifft_g pecos pecos_src ${TPL_LIBS} ${LAPACK_LIBS}
${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(pecos_ifft_g pecos_ifft_g)

# Complete these for each specific test.
add_executable(pecos_ifft_sd pecos_ifft_sd.cpp)
target_link_libraries(pecos_ifft_sd pecos pecos_src ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(pecos_ifft_sd pecos_ifft_sd)

# Complete these for each specific test.
add_executable(boost_test_dist boost_test_dist.cpp)
target_link_libraries(boost_test_dist pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(boost_test_dist boost_test_dist)

# Complete these for each specific test.
add_executable(boost_test_rng boost_test_rng.cpp)
target_link_libraries(boost_test_rng pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(boost_test_rng boost_test_rng)

# Complete these for each specific test.
add_executable(pecos_int_driver pecos_int_driver.cpp)
target_link_libraries(pecos_int_driver pecos pecos_src ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(pecos_int_driver pecos_int_driver)

# Complete these for each specific test.
add_executable(pecos_lhs_driver pecos_lhs_driver.cpp)
target_link_libraries(pecos_lhs_driver pecos pecos_src ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
add_test(pecos_lhs_driver pecos_lhs_driver)

