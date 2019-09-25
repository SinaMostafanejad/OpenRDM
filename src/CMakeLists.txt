add_executable(mcpdft main.cc
                      mcpdft.cc
                      energy.cc
                      functional.cc
                      utility.cc
		      libMem.cc)

target_include_directories(mcpdft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(mcpdft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mcpdft PUBLIC coverage_config)

install(TARGETS mcpdft
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR} 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBRID})
#install(TARGETS mcpdft
       #        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
       #        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
       #        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBRID})

