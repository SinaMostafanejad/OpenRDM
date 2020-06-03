Running the tests    {#runningtests}
=================

Since the current implementation of the <b>OpenRDM</b> package is managed via CMake,
the tests can be run using the "ctest" command. For seeing the results of the running
tests in the active terminal, add the -V option to the ctest command for a verbose output:

\code{.cmake}
ctest -V
\endcode

The list of available tests can be querried via adding the -N option to the ctest command:

\code{.cmake}
ctest -N
\endcode

The aforementioned command also prints the test number attributed to each test case that can
be individually selected and called via adding the -I option to the ctest command:

\code{.cmake}
ctest -I 2
\endcode

which in this case, runs the test \#2.
