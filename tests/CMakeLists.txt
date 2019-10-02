# list of test examples
list(APPEND TEST_SET h2_svwn_sto3g)
add_executable(mcpdft utility.cc
	              main.cc
)
target_link_libraries(mcpdft PUBLIC openrdm)
target_compile_definitions(mcpdft PUBLIC openrdm)
include_directories(${OPENRDM_SOURCE_DIR}/src)
if(WITH_LIBXC)
#  add_library(libxc_ext STATIC IMPORTED)
#  set_target_properties(libxc_ext PROPERTIES IMPORTED_LOCATION ${Libxc_LIBRARIES})
#  target_link_libraries(mcpdft PUBLIC libxc_ext)
   target_include_directories(mcpdft PRIVATE ${Libxc_INCLUDE_DIRS})
endif()
if(WITH_HDF5)
#  add_library(libxc_ext STATIC IMPORTED)
#  set_target_properties(libxc_ext PROPERTIES IMPORTED_LOCATION ${Libxc_LIBRARIES})
#  target_link_libraries(mcpdft PUBLIC libxc_ext)
   target_include_directories(mcpdft PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

foreach(test_id IN LISTS TEST_SET)
    add_test(NAME ${test_id}
	     COMMAND mcpdft "${test_id}"
	     WORKING_DIRECTORY ${OPENRDM_SOURCE_DIR}/tests
	    )
endforeach()
