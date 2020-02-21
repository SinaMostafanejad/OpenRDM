if(${WITH_LIBXC})
   find_package(Libxc 4.3.4 CONFIG QUIET)
   if(${Libxc_FOUND})
      message(STATUS "${Yellow}Libxc ${Libxc_VERSION} is found!${ColourReset}")
      message(STATUS "${Yellow}The list of Libxc libraries to be linked:${ColourReset} ${Libxc_LIBRARIES}")
      message(STATUS "${Yellow}The list of Libxc include directories:${ColourReset} ${Libxc_INCLUDE_DIRS}")
      #message(STATUS "${Yellow}The list of Libxc definitions:${ColourReset} ${Libxc_DEFINITIONS}")
   else()
      message(STATUS "${Cyan}Libxc is not found. Building from extrenal source...${ColourReset}")
      include(ExternalProject)
      ExternalProject_Add(Libxc
         GIT_REPOSITORY git@gitlab.com:libxc/libxc.git 
         GIT_TAG 31cf5534df849a8defc416a11e0b9938731d12e3
	 UPDATE_COMMAND ""
         DOWNLOAD_DIR ${OPENRDM_BINARY_DIR}/external/upstream/Libxc_download
         SOURCE_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/Libxc_src
         BINARY_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/Libxc_bin
	 INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/Libxc
         TMP_DIR      ${OPENRDM_BINARY_DIR}/external/upstream/Libxc_tmp
         STAMP_DIR    ${OPENRDM_BINARY_DIR}/external/upstream/Libxc_stamp
	 CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                    -DBUILD_SHARED_LIBS=OFF
		    -DBUILD_TESTING=OFF
                    -DBUILD_FPIC=ON
                    -DNAMESPACE_INSTALL_INCLUDEDIR=/
                    -DENABLE_GENERIC=OFF
                    -DENABLE_XHOST=ON
		    -DENABLE_FORTRAN=OFF
                    -DENABLE_FORTRAN03=OFF
		    -DDISABLE_VXC=ON
		    -DDISABLE_FXC=ON
                    -DDISABLE_KXC=ON
                    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                    -DCMAKE_INSTALL_LIBDIR=<INSTALL_DIR>/lib
                    -DCMAKE_INSTALL_BINDIR=<INSTALL_DIR>/bin
                    -DCMAKE_INSTALL_INCLUDEDIR=<INSTALL_DIR>/include
                    -DCMAKE_INSTALL_DATADIR=<INSTALL_DIR>/share
         BUILD_BYPRODUCTS=<INSTALL_DIR>/lib/libxc.a
         LOG_CONFIGURE YES
         LOG_BUILD YES
         LOG_INSTALL YES
         LOG_UPDATE YES
         LOG_OUTPUT_ON_FAILURE YES
         BUILD_ALWAYS YES
      )
      set(Libxc_DIR ${OPENRDM_STAGING_PREFIX}/Libxc/share/cmake/Libxc CACHE PATH "path to internally built LibxcConfig.cmake" FORCE)
      set(Libxc_INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/Libxc CACHE PATH "path to libxc root directory" FORCE)
      set(Libxc_INCLUDE_DIRS ${OPENRDM_STAGING_PREFIX}/Libxc/include CACHE PATH "path to libxc include directory" FORCE)
      set(Libxc_LIBRARIES    ${OPENRDM_STAGING_PREFIX}/Libxc/lib/libxc.a CACHE PATH "path to libxc library directory" FORCE)
      if(${Libxc_FOUND})
         message(STATUS "${Yellow}libxc ${Libxc_VERSION} is found!${ColourReset}")
         message(STATUS "${Yellow}The list of libxc libraries to be linked:${ColourReset} ${Libxc_LIBRARIES}")
         message(STATUS "${Yellow}The list of libxc include directories:${ColourReset} ${Libxc_INCLUDE_DIRS}")
         message(STATUS "${Yellow}The list of Libxc definitions:${ColourReset} ${Libxc_DEFINITIONS}")
         target_include_directories(openrdm INTERFACE
                                    $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/upstream/libxc_src/include>
                                    $<INSTALL_INTERFACE:${Libxc_INCLUDE_DIRS}>
         )
         target_link_libraries(openrdm PUBLIC ${Libxc_LIBRARIES})
      endif()
   endif()
endif()
