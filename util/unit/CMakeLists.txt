# CMakeLists.txt for util unit tests
#
include_directories(${Pecos_SOURCE_DIR}/util/src)
set(pecos_util_unit_tests
  linear_algebra_test.cpp
  linear_solvers_test.cpp
  math_tools_test.cpp
  )
add_executable(pecos_util_unit_tests 
  teuchos_unit_test_driver.cpp
  ${pecos_util_unit_tests}
  )

target_link_libraries(pecos_util_unit_tests pecos_util ${Teuchos_LIBRARIES})
add_test(pecos_util_unit_tests pecos_util_unit_tests)
