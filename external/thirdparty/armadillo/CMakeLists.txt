if(WITH_ARMADILLO)
#   include(FindArmadillo)
   if(ARMADILLO_FOUND)
      message(STATUS "${Yellow}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
      message(STATUS "${Yellow}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
      message(STATUS "${Yellow}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
   else()
      include(ExternalProject)
      ExternalProject_Add(armadillo
         GIT_REPOSITORY git@gitlab.com:conradsnicta/armadillo-code.git
         GIT_TAG 64c517324f996df27528ffdc519277c92b6ce0e7
         UPDATE_COMMAND ""
	 DOWNLOAD_DIR ${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_download
	 SOURCE_DIR   ${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_src
         BINARY_DIR   ${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_bin
	 INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/armadillo
	 TMP_DIR      ${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_tmp
         STAMP_DIR    ${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_stamp
	 CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
		    -DCMAKE_INSTALL_LIBDIR=<INSTALL_DIR>/lib
		    -DCMAKE_INSTALL_BINDIR=<INSTALL_DIR>/bin
		    -DCMAKE_INSTALL_INCLUDEDIR=<INSTALL_DIR>/include
		    -DCMAKE_INSTALL_DATADIR=<INSTALL_DIR>/share
         BUILD_ALWAYS YES
      )
      set(armadillo_DIR ${OPENRDM_STAGING_DIR}/armadillo/share/Armadillo/CMake CACHE PATH "path to internally built ArmadilloConfig.cmake" FORCE)
      set(armadillo_INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/armadillo CACHE PATH "path to armadillo root directory" FORCE)
      set(ARMADILLO_INCLUDE_DIRS ${OPENRDM_STAGING_PREFIX}/armadillo/include CACHE PATH "path to armadillo include directory" FORCE)
      set(ARMADILLO_LIBRARIES    ${OPENRDM_STAGING_PREFIX}/armadillo/lib/libarmadillo.so CACHE PATH "path to armadillo library directory" FORCE)
      message(STATUS "${Magenta}Armadillo ${ARMADILLO_VERSION_STRING} is found!${ColourReset}")
      message(STATUS "${Magenta}The list of Armadillo libraries to be linked:${ColourReset} ${ARMADILLO_LIBRARIES}")
      message(STATUS "${Magenta}The list of Armadillo include directories:${ColourReset} ${ARMADILLO_INCLUDE_DIRS}")
      #      if(ARMADILLO_FOUND)
      target_include_directories(openrdm PUBLIC 
                                    $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/thirdparty/armadillo_src/include>
                                    $<INSTALL_INTERFACE:${ARMADILLO_INCLUDE_DIRS}>
           			   )
			   target_link_libraries(openrdm PUBLIC ${ARMADILLO_LIBRARIES})
	 #      endif()
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
