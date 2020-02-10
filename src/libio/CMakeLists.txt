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
	    HDF5Client.cc)
target_sources(openrdm
               INTERFACE
	       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/${idx}>
	       $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
endforeach()

target_link_libraries(openrdm INTERFACE ${HDF5_LIBRARIES} dl)

target_include_directories(openrdm
                           INTERFACE
	                   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
	                   $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
