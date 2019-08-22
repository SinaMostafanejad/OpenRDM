cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(mcpdft LANGUAGES CXX C)

add_executable(mcpdft main.cc mcpdft.cc energy.cc)

set_property(TARGET mcpdft PROPERTY CXX_STANDARD 11)

#target_sources(mcpdft
#               PRIVATE
#               ${CMAKE_CURRENT_LIST_DIR}/mcpdft.cc
#               ${CMAKE_CURRENT_LIST_DIR}/mcpdft.h
#)
#
target_include_directories(mcpdft PRIVATE ${CMAKE_BINRAY_DIR})

