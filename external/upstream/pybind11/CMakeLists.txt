#find_package(pybind11 2.3.0 CONFIG QUIET)
#
#if(${pybind11_FOUND})
#    message(STATUS "${Cyan}Found pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
#    add_library(pybind11_external INTERFACE)  # dummy
#else()
#    if(${CMAKE_INSIST_FIND_PACKAGE_pybind11})
#        message(FATAL_ERROR "Suitable pybind11 could not be externally located as user insists")
#    endif()
#
#    include(ExternalProject)
#    message(STATUS "Suitable pybind11 could not be located, ${Magenta}Building pybind11${ColourReset} instead.")
#    ExternalProject_Add(pybind11_external
#        URL https://github.com/pybind/pybind11/archive/v2.3.0.tar.gz
#        UPDATE_COMMAND ""
#        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
#                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}  # ditto
#                   -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
#                   -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
#                   -DPYBIND11_TEST=OFF)
#
#    set(pybind11_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
#endif()

if(${WITH_PYBIND})
   add_dependencies(${PROJECT_NAME} pybind11)
   find_package(pybind11 CONFIG QUIET)
   if(${pybind11_FOUND})
      message(STATUS "${Yellow}pybind11 ${pybind11_VERSION} is found!${ColourReset}")
      message(STATUS "${Yellow}The list of pybind11 include directories:${ColourReset} ${pybind11_INCLUDE_DIRS}")
      target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
   else()
      message(STATUS "${Yellow}Building pybind11 from extrenal source...${ColourReset}")
      include(ExternalProject)
      ExternalProject_Add(pybind11
         GIT_REPOSITORY git@github.com:pybind/pybind11.git
         GIT_TAG 4f72ef846fe8453596230ac285eeaa0ce3278bb4
         UPDATE_COMMAND ""
         DOWNLOAD_DIR ${OPENRDM_BINARY_DIR}/external/upstream/pybind11_download
         SOURCE_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/pybind11_src
         BINARY_DIR   ${OPENRDM_BINARY_DIR}/external/upstream/pybind11_bin
         INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/pybind11
         TMP_DIR      ${OPENRDM_BINARY_DIR}/external/upstream/pybind11_tmp
         STAMP_DIR    ${OPENRDM_BINARY_DIR}/external/upstream/pybind11_stamp
         CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                    -DBUILD_SHARED_LIBS=OFF
                    -DPYBIND11_TEST=OFF
                    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                    -DCMAKE_INSTALL_LIBDIR=<INSTALL_DIR>/lib
                    -DCMAKE_INSTALL_BINDIR=<INSTALL_DIR>/bin
                    -DCMAKE_INSTALL_INCLUDEDIR=<INSTALL_DIR>/include
                    -DCMAKE_INSTALL_DATADIR=<INSTALL_DIR>/share
		    #BUILD_ALWAYS YES
      )
      set(pybind11_DIR ${OPENRDM_STAGING_DIR}/pybind11/share/cmake/pybind11 CACHE PATH "path to internally built pybind11Config.cmake" FORCE)
      set(pybind11_INSTALL_DIR  ${OPENRDM_STAGING_PREFIX}/pybind11 CACHE PATH "path to pybind11 root directory" FORCE)
      set(pybind11_INCLUDE_DIRS ${OPENRDM_STAGING_PREFIX}/pybind11/include CACHE PATH "path to pybind11 include directory" FORCE)
      if(${pybind11_FOUND})
         message(STATUS "${Yellow}The list of pybind11 libraries to be linked:${ColourReset} ${pybind11_LIBRARIES}")
         message(STATUS "${Yellow}The list of pybind11 include directories:${ColourReset} ${pybind11_INCLUDE_DIRS}")
         target_include_directories(${PROJECT_NAME} PUBLIC
                                    $<BUILD_INTERFACE:${OPENRDM_BINARY_DIR}/external/upstream/pybind11_src/include>
                                    $<INSTALL_INTERFACE:${pybind11_INCLUDE_DIRS}>
         )
      endif()
   endif()
endif()

