#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# CMakeLists.txt for util unit tests
#
include_directories(${Pecos_SOURCE_DIR}/util/src)

add_executable(pecos_linear_algebra_test  linear_algebra_test.cpp)
target_link_libraries(pecos_linear_algebra_test pecos_util ${Teuchos_LIBRARIES} Boost::boost)
add_test(pecos_linear_algebra_test pecos_linear_algebra_test)

add_executable(pecos_linear_solvers_test  linear_solvers_test.cpp)
target_link_libraries(pecos_linear_solvers_test pecos_util ${Teuchos_LIBRARIES} Boost::boost)
add_test(pecos_linear_solvers_test pecos_linear_solvers_test)

add_executable(pecos_math_tools_test  math_tools_test.cpp)
target_link_libraries(pecos_math_tools_test pecos_util ${Teuchos_LIBRARIES} Boost::boost)
add_test(pecos_math_tools_test pecos_math_tools_test)

add_executable(options_unit_tests  options.cpp)
target_link_libraries(options_unit_tests pecos_util ${Teuchos_LIBRARIES}
  Boost::boost)
add_test(options_unit_tests options_unit_tests)
