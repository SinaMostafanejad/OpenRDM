cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(mcpdft LANGUAGES CXX C)

add_subdirectory(src)

set_property(TARGET mcpdft PROPERTY CXX_STANDARD 11)
