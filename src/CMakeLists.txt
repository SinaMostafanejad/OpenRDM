set(sources_ "")
set(targets_ "${PROJECT_NAME}")
foreach(idx libmcpdft
	    libmem
            libio
            libfunc)
       list(APPEND sources_ "${idx}")
       add_subdirectory(${idx})
endforeach()
if(WITH_PSI4)
   list(APPEND sources_ libpsi4)
   list(APPEND targets_ mcpdft_iface)
   add_subdirectory(libpsi4)
endif()

set_target_properties(${PROJECT_NAME}
   PROPERTIES
      POSITION_INDEPENDENT_CODE 1
      ARCHIVE_OUTPUT_NAME "openrdm"
      DEBUG_POSTFIX "_sd"
      RELEASE_POSTFIX "_s"
      SKIP_BUILD_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_RPATH_USE_LINK_PATH ON
)
#-------------------------------------------------------------------
install(TARGETS
           ${targets_}
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
foreach(idx IN LISTS sources_)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${idx}/
               DESTINATION ${INSTALL_INCLUDEDIR}
               FILES_MATCHING
               PATTERN *.h
       )
endforeach()
#-------------------------------------------------------------------
unset(targets_)
unset(sources_)

