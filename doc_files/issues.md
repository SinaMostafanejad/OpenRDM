Known Issues     {#issues}
============

[TOC]

@section ninjaissue Ninja Compilers

Since Libxc and HDF5 are added as external projects and independent builds into the main <b>OpenRDM</b> project,
using <i>H\_LIBXC=ON</i> and <i>WITH\_HDF5=ON</i> options at the configure step generates
an "unknown build rule" error using Ninja generators. At the moment, the user has to adopt Unix Makefile's make
as the build generator.

@section hdf5issue HDF5 Version Inconsistency Error

Activating the <i>WITH\_HDF5=ON</i> CMake option triggers downloading HDF5 library from its repository
and installs it as an external project. If you already have the HDF5 library, header files or
header-only version installed on your system, cloning and installing a fresh copy as an external project
might cause a version inconsistency error.

In this situation, we recommend the user to uninstall/remove the older libraries and header files and
setup a fresh copy of the new version that matches that of the external HDF5 library cloned and 
installed by <b>OpenRDM</b>.

If the user cannot afford to refresh the local copy of HDF5 header-only library, adding the following
environment variable to the <i>.bashrc</i> file would supress the aforementioned error:

\code{.bash}
export HDF5_DISABLE_VERSION_CHECK=2
\endcode

We do not, however, recommend this way due to any possible instability issue or unpredicted behavior.
