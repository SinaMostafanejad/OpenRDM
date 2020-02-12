# list of test examples
list(APPEND TEST_SET h2_tsvwn_sto3g
	             h2_tpbe_sto3g)

add_executable(mcpdft main.cc)
target_link_libraries(mcpdft PUBLIC openrdm)
target_compile_definitions(mcpdft PUBLIC openrdm)

if(WITH_LIBXC)
   target_include_directories(mcpdft PRIVATE ${Libxc_INCLUDE_DIRS})
endif()

foreach(test_id IN LISTS TEST_SET)
    if("${test_id}" STREQUAL "h2_tsvwn_sto3g")
      set(functional_id "SVWN")
    else()
      set(functional_id "PBE")
    endif()
    add_test(NAME ${test_id}
             COMMAND mcpdft "${test_id}" ${functional_id}
             WORKING_DIRECTORY ${OPENRDM_SOURCE_DIR}/tests
	    )
endforeach()
