Known Issues     {#issues}
============

Since Libxc and HDF5 are added as external projects and independent builds into the main OpenRDM project,
using WITH_LIBXC=ON and WITH_HDF5=ON options at the configure step and within the configure script generates
an "unknown build rule" error for Ninja generators. At the moment, the user has to adopt Unix Makefile's make
as the build generator.
