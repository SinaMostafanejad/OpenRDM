To Do   {#todo}
=======

The following is a list of important tasks that need to be done:

+ MCPDFT:
   - Finish the implementation of the fully-translated density function.
+ HDF5 library:
   - Think about using creational patters instead of a inflexible interface
     class.
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

There are also long-term goals:

+ Numerical grid:
   - Need to write a code for numerical quadrature grid
     generation or consider interfacing with an open-
     source library and manipulating it.
