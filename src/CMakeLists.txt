
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}
                    ${PROJECT_SOURCE_DIR}/packages/teuchos/src
                    ${PROJECT_BINARY_DIR}/packages/teuchos/src
                    ${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)

## Pecos library sources.
SET(Pecos_src
  DataTransformation.cpp
  InverseTransformation.cpp
  FourierInverseTransformation.cpp
  KarhunenLoeveInverseTransformation.cpp
  SamplingInverseTransformation.cpp
  ProbabilityTransformation.cpp
  NatafTransformation.cpp
  BasisPolynomial.cpp
  LagrangeInterpPolynomial.cpp
  OrthogonalPolynomial.cpp
  ChebyshevOrthogPolynomial.cpp
  GenLaguerreOrthogPolynomial.cpp
  HermiteOrthogPolynomial.cpp
  JacobiOrthogPolynomial.cpp
  LaguerreOrthogPolynomial.cpp
  LegendreOrthogPolynomial.cpp
  NumericGenOrthogPolynomial.cpp
  SparseGridDriver.cpp
  LHSDriver.cpp
  pecos_stat_util.cpp
  pecos_data_types.hpp
  pecos_global_defs.hpp
  pecos_stat_util.hpp
  DataTransformation.hpp
  InverseTransformation.hpp
  FourierInverseTransformation.hpp
  KarhunenLoeveInverseTransformation.hpp
  SamplingInverseTransformation.hpp
  ProbabilityTransformation.hpp
  NatafTransformation.hpp
  BasisPolynomial.hpp
  LagrangeInterpPolynomial.hpp
  OrthogonalPolynomial.hpp
  ChebyshevOrthogPolynomial.hpp
  GenLaguerreOrthogPolynomial.hpp
  HermiteOrthogPolynomial.hpp
  JacobiOrthogPolynomial.hpp
  LaguerreOrthogPolynomial.hpp
  LegendreOrthogPolynomial.hpp
  NumericGenOrthogPolynomial.hpp
  SparseGridDriver.hpp
  LHSDriver.hpp
  BoostRNG_Monostate.hpp
)

IF(WITH_LHS)
  SET(Pecos_src ${Pecos_src} lhs_wrappers.f90)
  IF(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    SET_PROPERTY(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS -ffixed-form)
  ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    SET_PROPERTY(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS -fixed)
  ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
    SET_PROPERTY(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS -Mfixed)
  ENDIF()
ENDIF(WITH_LHS)

ADD_LIBRARY(pecos_src ${Pecos_src})
IF(WITH_LHS)
TARGET_LINK_LIBRARIES(pecos_src teuchos lhs ${LAPACK_LIBRARIES})
ELSE(WITH_LHS)
TARGET_LINK_LIBRARIES(pecos_src teuchos ${LAPACK_LIBRARIES})
ENDIF(WITH_LHS)

