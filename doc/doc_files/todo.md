To Do   {#todo}
=======

The following is a list of important tasks that need to be done:

+ HDF5 library:
   - Providing a simple example showing the application of
     HDF5 in a disk-based low-memory algorithm.
+ Libxc library:
   - Adopt Abstract Factory pattern to wrap the Libxc library
     and tune it towards MCPDFT in OpenRDM.
+ Documentation:
   - Notify the user of the absence of a CMake configure file
     that helps find_package() find the openrdm package. When provided,
     the expedient solution provided in "How to Use" section should be
     updated.
+ Interface:
   - OpenRDM needs have a (series) of class(es) for interfacing with QC softwares.
     For example, Psi4 provides intermediate files and PSIO class for dealing with this issue;
     Possibly, having a simple class would address the problem in Psi4. Other programs need to
     have a similar functionality to handle this problem. The next candidate in this direction
     would be PySCF.
   - The first PySCF example involves generating the integrals on the PySCF side and import OpenRDM
     functions to the python code through binding with pybind11. The second example shows how
     to store the integrals in HDF5 files using h5py and read them back into a C++ code for OpenRDM.
+ CMake:
   - Migrate to superbuild structure.
   - Manage the tests in a more elaborate way.
   - Provide installation and packaging procedure. The code is in need of a configure file for
     CMake find_package() to be able to be found by other projects.

There are also long-term goals:

+ Numerical grid:
   - Need to write a code for numerical quadrature grid
     generation or consider interfacing with an open-
     source library and manipulating it.
