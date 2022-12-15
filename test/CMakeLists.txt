#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

include_directories(${Pecos_SOURCE_DIR}/src
                    ${Pecos_SOURCE_DIR}/util/src
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

set(pecos_test_libs pecos_src lhs)

function(pecos_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${pecos_test_libs})
  add_test(${name} ${name})
endfunction()


# BMA: Historical rationale, prior to removing libpecos (empty.cxx). Leaving
# the explicit linkage to lhs to be safe, though removed libpecos. It seems
# just linking ${Pecos_LINK_LIBRARIES} would suffice.
#
# Ideally the test executables need to link to only "pecos" because it
# already links to its dependencies.  However, the lhs library defines
# at least one symbol (SSPEV) that is also provided by lapack.  The
# lhs library must come before lapack even though lhs does not link to
# lapack.  This ensures that the linker prefers symbols from lhs over
# conflicting symbols in lapack.

pecos_add_test(pecos_ifft_g)
pecos_add_test(pecos_ifft_sd)
set_tests_properties(pecos_ifft_sd PROPERTIES DEPENDS pecos_ifft_g)
pecos_add_test(boost_test_dist)
pecos_add_test(boost_test_rng)
pecos_add_test(pecos_int_driver)
pecos_add_test(pecos_gsg_driver)
pecos_add_test(pecos_lhs_driver)
pecos_add_test(pecos_pochhammer)
pecos_add_test(pecos_discrete_poly)
pecos_add_test(pecos_linear_solvers)
pecos_add_test(pecos_utils)

install(TARGETS pecos_gsg_driver DESTINATION bin)
