FIND_PACKAGE( PythonLibs )
FIND_PACKAGE( Numpy )

# Define the target for the swig runtime header
SET(SWIG_RUNTIME swigpyrun.h)
ADD_CUSTOM_COMMAND(
  OUTPUT ${SWIG_RUNTIME}
  COMMAND ${SWIG_EXECUTABLE} -python -external-runtime
  )

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${SWIG_RUNTIME})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${NumPy_INCLUDE_DIR})

include_directories( ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${DakotaSurrogates_SOURCE_DIR}/models/src
		    ${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
	            ${DakotaSurrogates_SOURCE_DIR}/util/src
		    ${Teuchos_INCLUDE_DIRS})

# add approximation module
SET(FOUNDATION_LIBRARIES pecos_util ${Teuchos_LIBRARIES})
#SET(SURROGATE_LIBRARIES pecos_wrapper models ${FOUNDATION_LIBRARIES} ${Teuchos_LIBRARIES})
SET(SURROGATE_LIBRARIES models ${FOUNDATION_LIBRARIES} )
set_source_files_properties( approximations.i PROPERTIES CPLUSPLUS ON )
set_source_files_properties( regression.i PROPERTIES CPLUSPLUS ON )
set_source_files_properties( math_tools.i PROPERTIES CPLUSPLUS ON )

set_source_files_properties( approximations.i PROPERTIES SWIG_FLAGS "")
set_source_files_properties( regression.i PROPERTIES SWIG_FLAGS "")
set_source_files_properties( math_tools.i PROPERTIES SWIG_FLAGS "")
swig_add_module( approximation python approximations.i )
swig_add_module( regression python regression.i )
swig_add_module( math_tools python math_tools.i )
swig_link_libraries( approximation ${SURROGATE_LIBRARIES} ${PYTHON_LIBRARIES})
swig_link_libraries( regression ${FOUNDATION_LIBRARIES} ${PYTHON_LIBRARIES})
swig_link_libraries( math_tools ${FOUNDATION_LIBRARIES} ${PYTHON_LIBRARIES})
