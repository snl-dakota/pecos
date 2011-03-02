# CMakeLists.txt for pecos "unit" tests

add_definitions("-DHAVE_CONFIG_H")

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${Teuchos_INCLUDE_DIRS})

if(HAVE_SPARSE_GRID)
  include_directories(${PROJECT_SOURCE_DIR}/packages/VPISparseGrid/src)
endif(HAVE_SPARSE_GRID)

# Ideally the test executables need to link to only "pecos" because it
# already links to its dependencies.  However, the lhs library defines
# at least one symbol (SSPEV) that is also provided by lapack.  The
# lhs library must come before lapack even though lhs does not link to
# lapack.  This ensures that the linker prefers symbols from lhs over
# conflicting symbols in lapack.
set(pecos_libs pecos pecos_src lhs)

# Complete these for each specific test.
add_executable(pecos_ifft_g pecos_ifft_g.cpp)
target_link_libraries(pecos_ifft_g ${pecos_libs})
add_test(pecos_ifft_g pecos_ifft_g)

# Complete these for each specific test.
add_executable(pecos_ifft_sd pecos_ifft_sd.cpp)
target_link_libraries(pecos_ifft_sd ${pecos_libs})
add_test(pecos_ifft_sd pecos_ifft_sd)

# Complete these for each specific test.
add_executable(boost_test_dist boost_test_dist.cpp)
target_link_libraries(boost_test_dist ${pecos_libs})
add_test(boost_test_dist boost_test_dist)

# Complete these for each specific test.
add_executable(boost_test_rng boost_test_rng.cpp)
target_link_libraries(boost_test_rng ${pecos_libs})
add_test(boost_test_rng boost_test_rng)

# Complete these for each specific test.
add_executable(pecos_int_driver pecos_int_driver.cpp)
target_link_libraries(pecos_int_driver ${pecos_libs})
add_test(pecos_int_driver pecos_int_driver)

# Complete these for each specific test.
add_executable(pecos_lhs_driver pecos_lhs_driver.cpp)
target_link_libraries(pecos_lhs_driver ${pecos_libs})
add_test(pecos_lhs_driver pecos_lhs_driver)

