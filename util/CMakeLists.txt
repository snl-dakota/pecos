# util is allowed to depend on Teuchos
include_directories(${Teuchos_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(unit)
