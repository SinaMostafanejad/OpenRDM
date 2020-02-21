if(${WITH_ARMADILLO})
   find_package(Armadillo 9.4 QUIET)
   add_dependencies(${PROJECT_NAME} Armadillo)
   if(${ARMADILLO_FOUND})
      message(STATUS "${Yellow}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
      message(STATUS "${Yellow}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
      message(STATUS "${Yellow}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
      target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
      target_link_libraries(${PROJECT_NAME} PUBLIC ${ARMADILLO_LIBRARIES})
   else()
      message(STATUS "${Yellow}Building Armadillo from external source...${ColourReset}")
      include(ExternalProject)
      ExternalProject_Add(Armadillo
         GIT_REPOSITORY git@gitlab.com:conradsnicta/armadillo-code.git
         GIT_TAG 64c517324f996df27528ffdc519277c92b6ce0e7
         UPDATE_COMMAND ""
	 DOWNLOAD_DIR ${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_download
	 SOURCE_DIR   ${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_src
         BINARY_DIR   ${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_bin
	 INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/Armadillo
	 TMP_DIR      ${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_tmp
         STAMP_DIR    ${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_stamp
	 CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
		    -DCMAKE_INSTALL_LIBDIR=<INSTALL_DIR>/lib
		    -DCMAKE_INSTALL_BINDIR=<INSTALL_DIR>/bin
		    -DCMAKE_INSTALL_INCLUDEDIR=<INSTALL_DIR>/include
		    -DCMAKE_INSTALL_DATADIR=<INSTALL_DIR>/share
         LOG_CONFIGURE YES
	 LOG_BUILD YES
         LOG_INSTALL YES
	 LOG_UPDATE YES
         LOG_OUTPUT_ON_FAILURE YES
         BUILD_ALWAYS YES
      )
      set(Armadillo_DIR ${OPENRDM_STAGING_PREFIX}/Armadillo/share/Armadillo/CMake CACHE PATH "path to internally built ArmadilloConfig.cmake" FORCE)
      set(Armadillo_INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/Armadillo CACHE PATH "path to armadillo root directory" FORCE)
      set(ARMADILLO_INCLUDE_DIRS ${OPENRDM_STAGING_PREFIX}/Armadillo/include CACHE PATH "path to armadillo include directory" FORCE)
      set(ARMADILLO_LIBRARIES    ${OPENRDM_STAGING_PREFIX}/Armadillo/lib/libarmadillo.so CACHE PATH "path to armadillo library directory" FORCE)
      if(${ARMADILLO_FOUND})
         message(STATUS "${Yellow}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
         message(STATUS "${Yellow}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
         message(STATUS "${Yellow}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
	 target_include_directories(${PROJECT_NAME} PUBLIC 
	                            $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_src/include>
                                    $<INSTALL_INTERFACE:${ARMADILLO_INCLUDE_DIRS}>
         )
	 target_link_libraries(${PROJECT_NAME} PUBLIC 
	                            $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/thirdparty/Armadillo_bin/libarmadillo.so>
				    $<INSTALL_INTERFACE:${ARMADILLO_LIBRARIES}>
         )
      endif()
   endif()
endif()

# I could use the more robust FetchContent module but
# it gives me less control over the installation directories etc.

#if(WITH_ARMADILLO)
#   include(FindArmadillo)
#   if(ARMADILLO_FOUND)
#      message(STATUS "${Cyan}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
#      message(STATUS "${Cyan}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
#      message(STATUS "${Cyan}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
#   else()
#      include(FetchContent)
#      FetchContent_Declare(armadillo
#         GIT_REPOSITORY git@gitlab.com:conradsnicta/armadillo-code.git
#         GIT_TAG 64c517324f996df27528ffdc519277c92b6ce0e7
#         SOURCE_DIR   ${CMAKE_CURRENT_BINARY_DIR}/armadillo_src
#         SOURCE_DIR   ${CMAKE_CURRENT_BINARY_DIR}/armadillo_src
#         BINARY_DIR   ${CMAKE_CURRENT_BINARY_DIR}/armadillo_bin
#         INSTALL_DIR  ${CMAKE_BINARY_DIR}/stage/armadillo
#         TMP_DIR      ${CMAKE_CURRENT_BINARY_DIR}/armadillo_tmp
#         STAMP_DIR    ${CMAKE_CURRENT_BINARY_DIR}/armadillo_stamp
#      )
#      FetchContent_GetProperties(armadillo)
#      if(NOT armadillo_POPULATED)
#	 FetchContent_Populate(armadillo)
#	 add_subdirectory(${armadillo_SOURCE_DIR} ${armadillo_BINARY_DIR})
#         message(STATUS "${Cyan}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
#         message(STATUS "${Cyan}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
#         message(STATUS "${Cyan}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
#      endif()
#   endif()
#endif()
