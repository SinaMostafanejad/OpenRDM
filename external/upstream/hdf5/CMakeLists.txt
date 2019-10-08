if(WITH_HDF5)
   include(FindHDF5)
   if(NOT HDF5_FOUND)
      message(STATUS "${Yellow}HDF5 ${HDF5_VERSION} is found!${ColourReset}")
      message(STATUS "${Yellow}The list of HDF5 libraries to be linked:${ColourReset} ${HDF5_LIBRARIES}")
      message(STATUS "${Yellow}The list of HDF5 include directories:${ColourReset} ${HDF5_INCLUDE_DIRS}")
      message(STATUS "${Yellow}The list of HDF5 definitions:${ColourReset} ${HDF5_DEFINITIONS}")
   else()
      include(ExternalProject)
      ExternalProject_Add(hdf5
         GIT_REPOSITORY https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git
         GIT_TAG 643d298823255a684862f4c6eb85ac5afde4b796 
         UPDATE_COMMAND ""
	 DOWNLOAD_DIR ${OPENRDM_BINARY_DIR}/external/upstream/hdf5_download
	 SOURCE_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/hdf5_src
         BINARY_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/hdf5_bin
	 INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/hdf5
	 TMP_DIR      ${OPENRDM_BINARY_DIR}/external/upstream/hdf5_tmp
         STAMP_DIR    ${OPENRDM_BINARY_DIR}/external/upstream/hdf5_stamp
	 CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
	            -DCMAKE_INSTALL_LIBDIR=<INSTALL_DIR>/lib
	            -DCMAKE_INSTALL_BINDIR=<INSTALL_DIR>/bin
	            -DCMAKE_INSTALL_INCLUDEDIR=<INSTALL_DIR>/include
	            -DCMAKE_INSTALL_DATADIR=<INSTALL_DIR>/share
         #BUILD_ALWAYS YES
      )
      set(HDF5_DIR ${OPENRDM_STAGING_DIR}/hdf5/share/cmake/hdf5 CACHE PATH "path to internally built HDF5Config.cmake" FORCE)
      set(HDF5_INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/hdf5 CACHE PATH "path to hdf5 root directory" FORCE)
      set(HDF5_INCLUDE_DIRS ${OPENRDM_STAGING_PREFIX}/hdf5/include CACHE PATH "path to hdf5 include directory" FORCE)
      set(HDF5_LIBRARIES    ${OPENRDM_STAGING_PREFIX}/hdf5/lib/libhdf5.a CACHE PATH "path to hdf5 library directory" FORCE)
      set(HDF5_HL_LIBRARIES ${OPENRDM_STAGING_PREFIX}/hdf5/lib/libhdf5_hl.a CACHE PATH "path to hdf5 library directory" FORCE)
      message(STATUS "${Cyan}HDF5 ${HDF5_VERSION} is found!${ColourReset}")
      message(STATUS "${Cyan}The list of HDF5 libraries to be linked:${ColourReset} ${HDF5_LIBRARIES}")
      message(STATUS "${Cyan}The list of HDF5 include directories:${ColourReset} ${HDF5_INCLUDE_DIRS}")
      if(HDF5_FOUND)
         target_include_directories(openrdm INTERFACE 
                                    $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/upstream/hdf5_src/src>
                                    $<INSTALL_INTERFACE:${HDF5_INCLUDE_DIRS}>
           			   )
         target_link_libraries(openrdm INTERFACE ${HDF5_LIBRARIES})
      endif()
   endif()
endif()
