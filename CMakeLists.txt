# The minimum required version for CMake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

cmake_policy(SET CMP0079 NEW)

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



# add compiler flags for code coverage
SET(CMAKE_CXX_FLAGS="-g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
SET(CMAKE_C_FLAGS="-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")

#SET(CMAKE_SHARED_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
SET(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")

target_link_libraries(mcpdft gcov)
