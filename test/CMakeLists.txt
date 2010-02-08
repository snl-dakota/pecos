SET(TPL_LIBS teuchos)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src
                    ${PROJECT_SOURCE_DIR}/packages/teuchos/src
                    ${PROJECT_BINARY_DIR}/packages/teuchos/src)

IF(WITH_LHS)
  SET(TPL_LIBS ${TPL_LIBS} lhs)
ENDIF(WITH_LHS)

IF(WITH_DFFTPACK)
  SET(TPL_LIBS ${TPL_LIBS} dfftpack)
ENDIF(WITH_DFFTPACK)

IF(WITH_FFTW)
  SET(TPL_LIBS ${TPL_LIBS} fftw3)
ENDIF(WITH_FFTW)

# Complete these for each specific test.
ADD_EXECUTABLE(pecos_ifft_g pecos_ifft_g.cpp)
TARGET_LINK_LIBRARIES(pecos_ifft_g pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS} pecos_src)
ADD_TEST(pecos_ifft_g pecos_ifft_g)

# Complete these for each specific test.
ADD_EXECUTABLE(pecos_ifft_sd pecos_ifft_sd.cpp)
TARGET_LINK_LIBRARIES(pecos_ifft_sd pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
ADD_TEST(pecos_ifft_sd pecos_ifft_sd)

# Complete these for each specific test.
ADD_EXECUTABLE(boost_test_dist boost_test_dist.cpp)
TARGET_LINK_LIBRARIES(boost_test_dist pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
ADD_TEST(boost_test_dist boost_test_dist)

# Complete these for each specific test.
ADD_EXECUTABLE(boost_test_rng boost_test_rng.cpp)
TARGET_LINK_LIBRARIES(boost_test_rng pecos ${TPL_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})
ADD_TEST(boost_test_rng boost_test_rng)

