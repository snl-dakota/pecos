SET(SOURCES PolynomialChaosExpansionWrap.cpp
            PolynomialChaosExpansion.cpp
            OrthogonalPolynomialBasis.cpp)
SET(HEADERS PolynomialChaosExpansionWrap.hpp
            PolynomialChaosExpansion.hpp
            OrthogonalPolynomialBasis.hpp)
add_library(pecos_wrapper ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
		    ${DakotaSurrogates_SOURCE_DIR}/models/src
            ${Pecos_SOURCE_DIR}/src
		    ${Teuchos_INCLUDE_DIRS})
if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
  list(APPEND Pecos_libs sparsegrid)
endif()


target_link_libraries(pecos_wrapper pecos_src pecos_util ${Teuchos_LIBRARIES})
