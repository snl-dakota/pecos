
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/teuchos/src
                    ${PROJECT_BINARY_DIR}/packages/teuchos/src)

## Pecos library sources.
SET(Pecos_src
  DataTransformation.cpp
  InverseTransformation.cpp
  FourierInverseTransformation.cpp
  KarhunenLoeveInverseTransformation.cpp
  SamplingInverseTransformation.cpp
  ProbabilityTransformation.cpp
  NatafTransformation.cpp
  pecos_stat_util.cpp
  LHSDriver.cpp
  BasisFunction.cpp
  DataTransformation.hpp
  InverseTransformation.hpp
  FourierInverseTransformation.hpp
  KarhunenLoeveInverseTransformation.hpp
  SamplingInverseTransformation.hpp
  ProbabilityTransformation.hpp
  NatafTransformation.hpp
  pecos_stat_util.hpp
  LHSDriver.hpp
  BasisFunction.hpp
)

IF(WITH_LHS)
  SET(Pecos_src ${Pecos_src} lhs_wrappers.f90)
  SET_PROPERTY(SOURCE lhs_wrappers.f90 PROPERTY COMPILE_FLAGS -ffixed-form)
ENDIF(WITH_LHS)

ADD_LIBRARY(pecos_src ${Pecos_src})
IF(WITH_LHS)
TARGET_LINK_LIBRARIES(pecos_src teuchos lhs ${LAPACK_LIBRARIES})
ELSE(WITH_LHS)
TARGET_LINK_LIBRARIES(pecos_src teuchos ${LAPACK_LIBRARIES})
ENDIF(WITH_LHS)

