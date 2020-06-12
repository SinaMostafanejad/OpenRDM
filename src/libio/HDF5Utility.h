#ifndef HDF5UTILITY_H
#define HDF5UTILITY_H

namespace mcpdft {

class HDF5Utility {
   public:
      /// Read the number of basis functions (AOs, MOs, NO, etc.)
      void read_nbfs(size_t &nao, size_t &nmo);

      /// Read active space variables (nactele, nactorb, ncore, etc.)
      void read_active_space_vars(size_t &nactele,
		                  size_t &nactorb,
				  size_t &ncore, 
				  size_t &nfrz, 
				  size_t &nocc, 
				  size_t &nvir);

      /// Read grids w, x, y, z
      void read_grids();
};

}

#endif // HDF5UTILITY_H
