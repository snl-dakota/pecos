SET(SOURCES PolynomialChaosExpansion.cpp)
SET(HEADERS PolynomialChaosExpansion.hpp)
add_library(pecos_wrapper ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
		    ${DakotaSurrogates_SOURCE_DIR}/models/src
		    ${Teuchos_INCLUDE_DIRS})

target_link_libraries(pecos_wrapper pecos_util ${Teuchos_LIBRARIES})
