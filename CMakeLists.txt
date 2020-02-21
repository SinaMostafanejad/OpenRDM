# The minimum required version for CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
# project name and its corresponding specifications
project(openrdm
       	LANGUAGES CXX C
	VERSION 0.2.0
)
set(${PROJECT_NAME}_URL "https://github.com/SinaMostafanejad/OpenRDM")
set(${PROJECT_NAME}_EMAIL "s_mostafanejad@fsu.edu")
set(${PROJECT_NAME}_LICENSE "BSD-3 Clause")
set(${PROJECT_NAME}_DESCRIPTION "An open-source library for reduced-density matrix-based analysis and computation")
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()
#===================================================================
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#===================================================================
add_library(openrdm STATIC "")
add_library(openrdm::openrdm ALIAS openrdm)
#===================================================================
# include modules from cmake directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(custom_color_messages)
include(ExternalProject)
#===================================================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "${Cyan}Build type set to${ColourReset} ${CMAKE_BUILD_TYPE}")
#-------------------------------------------------------------------
set(OPENRDM_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(OPENRDM_BINARY_DIR ${PROJECT_BINARY_DIR})
set(OPENRDM_STAGING_PREFIX ${OPENRDM_BINARY_DIR}/stage)
list(APPEND THIRDLIBLIST "")
list(APPEND EXTLIBLIST "")
#-------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "${OPENRDM_STAGING_PREFIX}" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "${Cyan}OpenRDM will be installed to${ColourReset} ${CMAKE_INSTALL_PREFIX}")
#===================================================================
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OPENRDM_STAGING_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OPENRDM_STAGING_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OPENRDM_STAGING_PREFIX}/${CMAKE_INSTALL_BINDIR})
#-------------------------------------------------------------------
# Offer the user the choice of overriding the installation directories
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake)
else()
  set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")
set(CMAKE_INSTALL_MESSAGE LAZY)
#-------------------------------------------------------------------
# Report to user
foreach(i LIB BIN INCLUDE CMAKE)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${i}DIR} path_ )
  message(STATUS "${Cyan}Installing ${i} components to${ColourReset} ${path_}")
  unset(path_)
endforeach()
#===================================================================
configure_file(
   "${OPENRDM_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.h.in"
   "${OPENRDM_SOURCE_DIR}/src/libmcpdft/${PROJECT_NAME}Config.h"
)
if(${WITH_HDF5})
   configure_file(
      "${OPENRDM_SOURCE_DIR}/cmake/hdf5Config.h.in"
      "${OPENRDM_SOURCE_DIR}/src/libio/hdf5Config.h"
   )
endif()
#===================================================================
option(CODE_COVERAGE "${Cyan}Enable coverage reporting${ColourReset}" OFF)
option(BUILD_TESTING "${Cyan}Build the testing tree${ColourReset}" ON)
option(WITH_ARMADILLO "${cyan}Build with Armadillo linear algebra library${ColourReset}" ON)
option(WITH_LIBXC "${cyan}Build with Libxc density functional library${ColourReset}" OFF)
option(WITH_HDF5 "${cyan}Build with HDF5${ColourReset}" ON)
option(WITH_PYBIND "${cyan}Build with pybind11${ColourReset}" OFF)
option(WITH_OPENMP "${cyan}Build with OpenMP${ColourReset}" OFF)
option(WITH_PSI4 "${cyan}Enable Psi4 interface${ColourReset}" OFF)
if(${WITH_HDF5})
	option(WITH_CHUNKING "${cyan}compile libio with chunking enabled in HDF5${ColourReset}" ON)
	if(${WITH_CHUNKING})
           option(CHUNK_DIM0 "${cyan}HDF5 chunk size (dim 0)${ColourReset}" "1")
           option(CHUNK_DIM1 "${cyan}HDF5 chunk size (dim 1)${ColourReset}" "1")
	endif()
endif()
#===================================================================
if(${WITH_PSI4})
   find_package(PythonInterp 3.0 REQUIRED)
   find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EXACT REQUIRED)
   find_package(psi4 1.1 CONFIG REQUIRED)
   if(psi4_FOUND)
      get_property(loc_ TARGET psi4::core PROPERTY LOCATION)
      message(STATUS "${Yellow}Psi4 (${psi4_VERSION}) is found in ${loc_} ${ColourReset}")
      unset(loc_)
   endif()
endif()
#-------------------------------------------------------------------
if(${WITH_OPENMP})
   find_package(OpenMP REQUIRED QUIET)
   if(${OpenMP_FOUND})
      target_link_libraries(${PROJECT_NAME}
                            PRIVATE
		            OpenMP::OpenMP_CXX
      )
   else()
      message(FATAL_ERROR "OpenMP is not found!" )
   endif()
endif()
#===================================================================
# Code Coverage Configuration
#add_library(coverage_config INTERFACE)
#
#if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#  # Add required flags (GCC & LLVM/Clang)
#  target_compile_options(coverage_config INTERFACE
#    -O0        # no optimization
#    -g         # generate debug info
#    --coverage # sets all required flags
#  )
#  set(CMAKE_CXX_FLAGS="-g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
#  set(CMAKE_C_FLAGS="-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
#  set(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
#  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
#  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
#    target_link_options(coverage_config INTERFACE --coverage)
#  else()
#    target_link_libraries(coverage_config INTERFACE --coverage)
#  endif()
#endif()
#===================================================================
add_subdirectory(external)
add_subdirectory(src)
#===================================================================
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
include(CPack)
#===================================================================
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
# enabling ctest command for running test in the build tree
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
  if(${WITH_LIBXC})
     add_dependencies(mcpdft libxc)
  endif()
  if(${WITH_HDF5})
     add_dependencies(mcpdft hdf5)
  endif()
endif()
#===================================================================
