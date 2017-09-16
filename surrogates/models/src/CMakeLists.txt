SET(SOURCES Function.cpp CppFunction.cpp
	    Approximation.cpp	
	    PolyApproximation.cpp Monomial.cpp
            Variables.cpp BoundedVariables.cpp
	    VariableTransformation.cpp AffineVariableTransformation.cpp
	    RegressionBuilder.cpp)
SET(HEADERS Function.hpp CppFunction.hpp
	    Approximation.hpp
	    PolyApproximation.hpp Monomial.hpp
            Variables.hpp BoundedVariables.hpp
	    VariableTransformation.hpp AffineVariableTransformation.hpp
	    RegressionBuilder.hpp)
add_library(models ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/models/src
	            ${Teuchos_INCLUDE_DIRS})

target_link_libraries(models pecos_util ${Teuchos_LIBRARIES})
