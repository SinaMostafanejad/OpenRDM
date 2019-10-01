# list of test examples
list(APPEND TEST_SET h2_svwn_sto3g)
add_executable(mcpdft utility.cc
	              main.cc
)		

target_link_libraries(mcpdft PUBLIC openrdm)
target_compile_definitions(mcpdft PUBLIC openrdm)
include_directories(${CMAKE_SOURCE_DIR}/src)

foreach(test_id IN LISTS TEST_SET)
    add_test(NAME ${test_id}
	     COMMAND mcpdft "${test_id}"
	     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests
	    )
	    
endforeach()
