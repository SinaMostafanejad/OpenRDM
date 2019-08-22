cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(mcpdft LANGUAGES CXX C)

add_executable(mcpdft main.cc mcpdft.cc energy.cc)

set_property(TARGET mcpdft PROPERTY CXX_STANDARD 11)
