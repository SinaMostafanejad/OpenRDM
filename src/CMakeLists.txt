add_library(openrdm INTERFACE)

foreach(idx libmcpdft
	    libmem
            libio
            libfunc)
       add_subdirectory(${idx})
endforeach()

target_include_directories(openrdm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libmcpdft>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libmem>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libio>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libfunc>
  $<INSTALL_INTERFACE:${CMAKE_STAGING_PREFIX}/openrdm/include
)

target_link_libraries(openrdm INTERFACE coverage_config)

install(TARGETS openrdm
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBRID})

