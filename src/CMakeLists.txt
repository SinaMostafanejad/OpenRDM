add_executable(mcpdft
               main.cc
               mcpdft.cc
               energy.cc
               functional.cc
               utility.cc)

       #install(TARGETS mcpdft
       #        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
       #        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
       #        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBRID})

install(TARGETS mcpdft
	RUNTIME DESTINATION ${CMAKE_SOURCE_DIR} 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBRID})
