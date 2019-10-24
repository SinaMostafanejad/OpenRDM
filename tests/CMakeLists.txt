# list of test examples
list(APPEND TEST_SET h2_tsvwn_sto3g
	             h2_tpbe_sto3g)

add_executable(mcpdft main.cc
	              utility.cc
)
target_link_libraries(mcpdft PUBLIC openrdm)
target_compile_definitions(mcpdft PUBLIC openrdm)
#target_include_directories(mcpdft PUBLIC ${OPENRDM_SOURCE_DIR}/src/libmcpdft)

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
#target_include_directories(mcpdft PRIVATE ${HDF5_INCLUDE_DIRS} dl)
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
