How to Install    {#installation}
==============

The configuration, generation, build, installation, packaging and deployment phases of the
<b>OpenRDM</b> project are managed by the CMake build system generator. CMake provides a
great amount of control over the fine-tuning of various (optional) capabilities available in
the project at the target level. In this way, users can benefit from extensibility and
portability of the code yet easily be able to recast it towards their needs.

For the sake of user's convenience, we have provided a <i>configure</i> script
in the source directory reflecting the aforementioned phases to highlight the 
available options provided with the users in order to add/remove various components
and dependencies to the OpenRDM project.

In case any of the mandatory [dependencies](https://sinamostafanejad.github.io/OpenRDM/dependencies.html)
such as armadillo linear algebra package could not be found by <b>OpenRDM</b>, it fetches it from the
corresponding host website and builds it as an external project on the host disk.

By default, <b>OpenRDM</b> builds and installs itself internally within the <i>build/</i>
folder at its source root directory. Feel free to install it at any other location on your 
system by adding

\code{.cmake}
...
-DCMAKE_INSTALL_PREFIX="/your/new/installation/address"
...
\endcode

to the <i>configure</i> script in the root directory of the <b>OpenRDM</b> project.
