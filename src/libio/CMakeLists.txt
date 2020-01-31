target_sources(openrdm
               INTERFACE
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Contiguous.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Contiguous.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Contiguous.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Contiguous.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Compact.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Compact.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Compact.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Compact.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Chunked.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Read_Chunked.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Chunked.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5_Write_Chunked.h
	       ${CMAKE_CURRENT_LIST_DIR}/IWrite.h
	       ${CMAKE_CURRENT_LIST_DIR}/IRead.h
	       ${CMAKE_CURRENT_LIST_DIR}/IOFactory.h
	       #	       ${CMAKE_CURRENT_LIST_DIR}/IOFactory.cc
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5ContiguousFactory.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5CompactFactory.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5ChunkedFactory.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5Client.h
	       ${CMAKE_CURRENT_LIST_DIR}/HDF5Client.cc
)

target_include_directories(openrdm INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(openrdm INTERFACE ${HDF5_LIBRARIES} dl)
