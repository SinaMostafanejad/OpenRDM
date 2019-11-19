How to Use    {#howtouse}
===========

The present implementation of the <b>OpenRDM</b> package is a header-only library.
It can be easily incorporated into a project using CMake. Add the following commands
to the CMakeLists.txt file in the root directory of your software:

\code{.cmake}
add_executable(myTarget myCode.cc)
target_link_libraries(myTarget PUBLIC openrdm)
target_compile_definitions(myTarget PUBLIC openrdm)
\endcode
