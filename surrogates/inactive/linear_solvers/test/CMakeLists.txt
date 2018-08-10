# CMakeLists.txt for linear_solvers "unit" tests

include_directories(${DakotaSurrogates_SOURCE_DIR}/math_utils/src
                    ${DakotaSurrogates_SOURCE_DIR}/linear_solvers/src
                    ${Boost_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS})

remove_definitions("-DHAVE_CONFIG_H")

option(BOOST_DISABLE_ASSERTS "Toggle BoostAssert, default is AssertDisabled" ON)
if (BOOST_DISABLE_ASSERTS)
  add_definitions("-DBOOST_DISABLE_ASSERTS")
endif(BOOST_DISABLE_ASSERTS)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "XL" OR
    ${CMAKE_CXX_COMPILER_ID} MATCHES "VisualAge")
  add_definitions("-DBOOST_MPL_CFG_ASSERT_USE_RELATION_NAMES")
endif()


if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
endif(HAVE_SPARSE_GRID)

set(linear_solver_tests
  linear_solvers.cpp
)

add_executable(linear_solver_tests teuchos_unit_test_driver.cpp ${linear_solver_tests} )
target_link_libraries(linear_solver_tests linear_solvers math_utils ${Teuchos_LIBRARIES})
add_test(linear_solver_tests linear_solver_tests)
