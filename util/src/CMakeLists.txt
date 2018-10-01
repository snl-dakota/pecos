# Pecos utilities that should depend only on STL and Teuchos

set(UTIL_SOURCES
  least_angle_regression.cpp
  linear_algebra.cpp
  linear_solvers.cpp
  math_tools.cpp
  sparse_matrices.cpp
  CrossValidationIterator.cpp
  LinearSystemCrossValidationIterator.cpp
  LSQCrossValidationIterator.cpp
  CrossValidatedSolver.cpp
)

set(UTIL_HEADERS
  least_angle_regression.hpp
  linear_algebra.hpp
  LinearSystemSolver.hpp
  linear_solvers.hpp
  math_tools.hpp
  orthogonal_matching_pursuit.hpp
  sparse_matrices.hpp
  CrossValidationIterator.hpp
  LinearSystemCrossValidationIterator.hpp
  LSQCrossValidationIterator.hpp
  OMPSolver.hpp
  LARSolver.hpp	
  LSQSolver.hpp
  EqConstrainedLSQSolver.hpp
  CrossValidatedSolver.hpp
  )

add_library(pecos_util ${UTIL_SOURCES} ${UTIL_HEADERS})
#SET(CMAKE_INCLUDE_CURRENT_DIR ON)
#target_link_libraries(pecos_util ${Teuchos_LIBRARIES})
install(TARGETS pecos_util EXPORT ${ExportTarget} DESTINATION lib)
