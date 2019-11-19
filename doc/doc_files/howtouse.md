How to Use    {#howtouse}
===========

The present implementation of the <b>OpenRDM</b> package is a header-only library.
It can be easily incorporated into a project using CMake. To do so, add the following
commands to the CMakeLists.txt file in the root directory of your software:

\code{.cmake}
\# A short block of the CMake code needed to link the OpenRDM library \#
add_executable(myTarget myCode.cc)
target_link_libraries(myTarget PUBLIC openrdm)
target_compile_definitions(myTarget PUBLIC openrdm)
\endcode

The general structure of a driver code adopting some of the functionalities provided by various
classes of the <b>OpenRDM</b> library can be found in the "tests" folder of its root directory.
More complicated test samples are comming soon!
