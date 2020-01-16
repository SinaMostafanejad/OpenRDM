# The minimum required version for CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
# project name and its corresponding specifications
project(openrdm LANGUAGES CXX C)
set(mydft_URL "https://github.com/SinaMostafanejad/OpenRDM")
set(mydft_EMAIL "s_mostafanejad@fsu.edu")
set(mydft_LICENSE "BSD-3 Clause")
set(mydfft_DESCRIPTION "An open-source library for reduced-density matrix-based analysis and computation")
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()
#===================================================================
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#===================================================================
# include modules from cmake directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)
include(custom_color_messages)
include(ExternalProject)
#===================================================================
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "${Red}OpenRDM install:${ColourReset} ${CMAKE_INSTALL_PREFIX}")
set(OPENRDM_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(OPENRDM_BINARY_DIR ${PROJECT_BINARY_DIR})
set(OPENRDM_STAGING_PREFIX ${OPENRDM_BINARY_DIR}/stage)
list(APPEND THIRDLIBLIST "")
list(APPEND EXTLIBLIST "")
#===================================================================
configure_file(
   "${OPENRDM_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.h.in"
   "${OPENRDM_SOURCE_DIR}/src/libmcpdft/${PROJECT_NAME}Config.h"
)
configure_file(
   "${OPENRDM_SOURCE_DIR}/cmake/hdf5Config.h.in"
   "${OPENRDM_SOURCE_DIR}/src/libio/hdf5Config.h"
)
#===================================================================
option(CODE_COVERAGE "${Cyan}Enable coverage reporting${ColourReset}" OFF)
option(BUILD_TESTING "${Cyan}Build the testing tree${ColourReset}" ON)
option(WITH_ARMADILLO "${cyan}Build with Armadillo linear algebra library${ColourReset}" ON)
option(WITH_LIBXC "${cyan}Build with Libxc density functional library${ColourReset}" OFF)
option(WITH_HDF5 "${cyan}Build with HDF5${ColourReset}" ON)
option(WITH_OPENMP "${cyan}Build with HDF5${ColourReset}" OFF)
if(${WITH_HDF5})
   option(ENABLE_CHUNKING "${cyan}compile libio with chunking enabled in HDF5${ColourReset}" ON)
endif()
option(CHUNK_DIM0 "${cyan}HDF5 chunk size (dim 0)${ColourReset}" "1")
option(CHUNK_DIM1 "${cyan}HDF5 chunk size (dim 1)${ColourReset}" "1")
#===================================================================
# Code Coverage Configuration
add_library(coverage_config INTERFACE)

if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  set(CMAKE_CXX_FLAGS="-g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS="-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()
#===================================================================
add_subdirectory(src)
add_subdirectory(external)
if(WITH_ARMADILLO)
   add_dependencies(openrdm armadillo)
endif()
if(WITH_LIBXC)
   add_dependencies(openrdm libxc)
endif()
if(WITH_HDF5)
   add_dependencies(openrdm hdf5)
endif()
#===================================================================
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
# enabling ctest command for running test in the build tree
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
  if(WITH_LIBXC)
     add_dependencies(mcpdft libxc)
  endif()
  if(WITH_HDF5)
     add_dependencies(mcpdft hdf5)
  endif()
endif()
#===================================================================
