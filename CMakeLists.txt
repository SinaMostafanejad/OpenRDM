# The minimum required version for CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# project name and its corresponding specifications
project(mcpdft LANGUAGES CXX C)

# enabling ctest command for running test in the build tree
enable_testing()

# include modules from cmake directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)

# add a new subdirectory to the build tree
add_subdirectory(src)
add_subdirectory(tests)

set_property(TARGET mcpdft PROPERTY CXX_STANDARD 11)
