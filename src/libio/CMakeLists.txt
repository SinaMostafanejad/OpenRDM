foreach(idx HDF5_Read_Contiguous.cc
	    HDF5_Read_Contiguous.h
	    HDF5_Write_Contiguous.cc
	    HDF5_Write_Contiguous.h
	    HDF5_Read_Compact.cc
	    HDF5_Read_Compact.h
	    HDF5_Write_Compact.cc
	    HDF5_Write_Compact.h
	    HDF5_Read_Chunked.cc
	    HDF5_Read_Chunked.h
	    HDF5_Write_Chunked.cc
	    HDF5_Write_Chunked.h
	    IWrite.h
	    IRead.h
	    IOFactory.h
	    HDF5ContiguousFactory.h
	    HDF5CompactFactory.h
	    HDF5ChunkedFactory.h
	    HDF5Client.h
	    HDF5Client.cc
	    HDF5ErrorManager.cc
	    HDF5ErrorManager.h
	    HDF5Utility.cc
	    HDF5_Read_Sparse.cc
	    HDF5Utility.h)
target_sources(openrdm
               PRIVATE
	       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
	       $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
endforeach()

target_include_directories(openrdm
                           PUBLIC
	                   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
	                   $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

target_include_directories(openrdm PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(openrdm PUBLIC ${ARMADILLO_LIBRARIES})

target_link_libraries(openrdm PUBLIC ${HDF5_LIBRARIES} ${CMAKE_DL_LIBS})
