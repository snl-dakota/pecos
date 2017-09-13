SET(SOURCES least_angle_regression.cpp compressed_sensing.cpp)
SET(HEADERS LinearSystemSolver.hpp least_angle_regression.hpp orthogonal_matching_pursuit.hpp compressed_sensing.hpp)
add_library(linear_solvers ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/math_utils/src
                    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(math_utils)
