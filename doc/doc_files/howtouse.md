How to Use    {#howtouse}
===========

The present implementation of the <b>OpenRDM</b> project is packaged as a static library.
It can be easily consumed by a host project via the CMake build system generator. 
To do so, after the installation of OpenRDM in, let's say, <i>/home/OpenRDM</i> directory,
add the following commands to the <i>CMakeLists.txt</i> file in the root directory of the host 
(target) software:

\code{.cmake}
# Set the openrdm_DIR variable to the OpenRDM's "openrdmConfig.cmake" config file address
set(openrdm_DIR "home/OpenRDM/build/stage/share/cmake/openrdm")

# Find the OpenRDM library
find_package(openrdm CONFIG REQUIRED)

# Add the following to see whether OpenRDM is found by find_package()
if(openrdm_FOUND)
   message(STATUS "OpenRDM (${openrdm_VERSION}) is found in: ${openrdm_DIR}")
   # Link the OpenRDM to the target project
   target_link_libraries(${PROJECT_NAME} PRIVATE openrdm::openrdm)
endif()
\endcode

Let us know on [Slack](https://openrdm.slack.com/join/shared_invite/enQtOTM2MDg2MzUxNjIyLWNlMzFlOWFhYTVlZGQ3ZGYxNWY3NTk4ZjRhYzM3MTU5MWZhN2VhY2Y5NzBiNjVjYzU1YWJkZDc2ODdhYTM4Yjg) or open an issue on the [main page](https://github.com/SinaMostafanejad/OpenRDM) of the <b>OpenRDM</b> repository
should you have any problems or need any assistance regarding to the <b>OpenRDM</b> project.
