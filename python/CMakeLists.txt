configure_file("setup.py.cmake.in" "setup.py" @ONLY)

set(pydakota_stamp "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

add_custom_command(
  OUTPUT ${pydakota_stamp}
  COMMAND ${PYTHON_EXECUTABLE} setup.py build
  COMMAND ${CMAKE_COMMAND} -E touch ${pydakota_stamp}
##  DEPENDS ${DEPS}
  )

add_custom_target(pydakota DEPENDS ${pydakota_stamp})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
  "build;swigpyrun.h")

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} setup.py install)")
