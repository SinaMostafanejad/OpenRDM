target_sources(openrdm
               INTERFACE
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write.h
	       ${CMAKE_CURRENT_LIST_DIR}/IWrite.h
	       ${CMAKE_CURRENT_LIST_DIR}/IRead.h
	       ${CMAKE_CURRENT_LIST_DIR}/IOFactory.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5Factory.h
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(openrdm INTERFACE ${HDF5_LIBRARIES} dl)
