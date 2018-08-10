# CMakeLists.txt for math_utils unit tests
#
include_directories(${DakotaSurrogates_SOURCE_DIR}/math_utils/src)
set(surrogate_teuchos_unit_tests
  linear_algebra_test.cpp
  math_tools_test.cpp
  random_vars.cpp
  )
add_executable(surrogate_teuchos_unit_tests 
  teuchos_unit_test_driver.cpp
  ${surrogate_teuchos_unit_tests}
  )

target_link_libraries(surrogate_teuchos_unit_tests math_utils ${Teuchos_LIBRARIES})
add_test(surrogate_teuchos_unit_tests surrogate_teuchos_unit_tests)
