configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Multivariate Function Approximation\ With Dakota.ipynb.in" "${CMAKE_CURRENT_BINARY_DIR}/Multivariate Function Approximation\ With Dakota.ipynb" COPY_ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_notebooks.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/test_notebooks.py COPY_ONLY)
