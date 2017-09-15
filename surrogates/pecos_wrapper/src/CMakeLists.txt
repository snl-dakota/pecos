SET(SOURCES PolynomialChaosExpansionWrap.cpp)
SET(HEADERS PolynomialChaosExpansionWrap.hpp)
add_library(pecos_wrapper ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
		    ${DakotaSurrogates_SOURCE_DIR}/models/src
            ${Pecos_SOURCE_DIR}/src
		    ${Teuchos_INCLUDE_DIRS})

target_link_libraries(pecos_wrapper pecos_src pecos_util ${Teuchos_LIBRARIES})
