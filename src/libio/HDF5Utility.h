#ifndef HDF5UTILITY_H
#define HDF5UTILITY_H

namespace mcpdft {

class HDF5Utility {
   public:
      /// Read the number of basis functions (AOs, MOs, NO, etc.)
      void read_nbfs(int &nao, int &nmo);

      /// Read grids w, x, y, z
      void read_grids();
};

}

#endif // HDF5UTILITY_H
