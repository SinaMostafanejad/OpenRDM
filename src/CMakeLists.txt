foreach(idx libmcpdft
	    libmem
            libio
            libfunc)
       add_subdirectory(${idx})
endforeach()

#target_link_libraries(openrdm INTERFACE coverage_config)
#-------------------------------------------------------------------
install(TARGETS
           openrdm
	EXPORT
           ${PROJECT_NAME}_Targets
	ARCHIVE
	   DESTINATION ${INSTALL_LIBDIR}
	   COMPONENT lib
        LIBRARY
           DESTINATION ${INSTALL_LIBDIR}
	   COMPONENT lib
	PUBLIC_HEADER
	DESTINATION ${INSTALL_INCLUDEDIR}
	   COMPONENT dev
#        INCLUDES ${CMAKE_INSTALL_INCLUDEDIR}
)
write_basic_package_version_file(
	"${OPENRDM_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
   "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
   INSTALL_DESTINATION ${INSTALL_INCLUDEDIR}
)
install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE "${PROJECT_NAME}::"
        DESTINATION ${INSTALL_CMAKEDIR}
)
#-------------------------------------------------------------------
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${INSTALL_CMAKEDIR}
)
foreach(idx libmcpdft
	    libmem
            libio
            libfunc)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${idx}/
               DESTINATION ${INSTALL_INCLUDEDIR}
               FILES_MATCHING
               PATTERN *.h
	       PATTERN *.cc
       )
endforeach()
