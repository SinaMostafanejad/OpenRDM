@page todo To Do 

The following is a list of important tasks that need to be done:

+ MCPDFT:
   - Finish the implementation of the fully-translated density function.
+ HDF5 library:
   - Interfacing the HDF5 with OpenRDM within the diskRW.cc implementation.
   - Providing a simple example showing the application of
     HDF5 in a disk-based low-memory algorithm.
+ Libxc library:
   - Adopt some design patterns to wrap the Libxc library
     and tune it towards MCPDFT in OpenRDM.
+ Documentation:
   - Provide a concise and clear explanation of MCPDFT and
     a more detailed documentation for the member variables
     and functions.
   - Consider migration from Doxygen to Sphinx or keeping
     both which might be redundant.

There are also long-term goals:

+ Numerical grid:
   - Need to write a code for numerical quadrature grid
     generation or consider interfacing with an open-
     source library and manipulating it.
