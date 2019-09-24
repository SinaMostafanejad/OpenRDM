# The minimum required version for CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# project name and its corresponding specifications
project(mcpdft LANGUAGES CXX C)
#===================================================================
# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)
#===================================================================
# include modules from cmake directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)
#===================================================================
option(CODE_COVERAGE "${Cyan}Enable coverage reporting${ColourReset}" OFF)
option (BUILD_TESTING "${Cyan}Build the testing tree.${ColourReset}" ON)
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
  SET(CMAKE_CXX_FLAGS="-g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
  SET(CMAKE_C_FLAGS="-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
  #SET(CMAKE_SHARED_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
  SET(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()

# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
# enabling ctest command for running test in the build tree
  enable_testing()
  add_subdirectory(tests)
endif()
#===================================================================
# add a new subdirectory to the build tree
add_subdirectory(src)
set_property(TARGET mcpdft PROPERTY CXX_STANDARD 11)
