Installation    {#installation}
============

The configuration, generation, build, installation, packaging and deployment phases of the
<b>OpenRDM</b> project is managed by the CMake build system generator. Cmake provides a
great amount of control over fine-tuning of various (optional) capabilities available in
the project at the target level. In this way, users can benefit from extensibility and
portability of the code yet easily be able to recast it towards their needs.

For the sake of user's convenience, we have provided a <i>configure</i> script
in the source directory to highlight the available options for adding/removing
various components and dependencies to the OpenRDM project.

In case any of the mandatory [dependencies](https://sinamostafanejad.github.io/OpenRDM/dependencies.html)
such as armadillo linear algebra package cannot be found by <b>OpenRDM</b>, it will fetch it from its
host website and builds it as an external project on your disk.

By default, <b>OpenRDM</b> builds and install itself in the <i>build/</i> folder at its source root
directory. Feel free to install it at any other location by adding 

\code{.cmake}
...
-DCMAKE_INSTALL_PREFIX="/your/new/installation/address"
...
\endcode

to the <i>configure</i> script in the root directory of the <b>OpenRDM</b> project.
