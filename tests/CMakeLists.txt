# list of test examples
set(TEST_SET h2_svwn_sto3g)

foreach(TEST_SET test_id)
    add_test(NAME h2-svwn-sto3g
             COMMAND mcpdft "${test_id}"
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
