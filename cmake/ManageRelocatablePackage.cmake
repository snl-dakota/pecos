# Manage a package that may appear in the local packages/ directory,
# or may have been previously defined with PACKAGE_DIR by a parent
#   package: package name in capitals
#   dirname: directory containing package 
function(ManageRelocatablePackage package dirname)
  if(HAVE_${package})
    message(STATUS
      "Finding relocatable package ${package} with directory ${dirname}")
    if(NOT ${package}_DIR)
      if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dirname}")
        set(${package}_DIR "${CMAKE_CURRENT_BINARY_DIR}/${dirname}")
	message(STATUS 
	  "Setting ${package}_DIR = ${CMAKE_CURRENT_BINARY_DIR}/${dirname}")
	add_subdirectory(${dirname})
      else()
	message(SEND_ERROR
	  "Directory ${CMAKE_CURRENT_SOURCE_DIR}/${dirname} missing")
      endif()
    endif()
  endif()
endfunction()
