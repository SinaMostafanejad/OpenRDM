#ifndef HDF5_ERROR_MANAGER_H
#define HDF5_ERROR_MANAGER_H

#include <string>
#include "hdf5.h"

namespace mcpdft {

class HDF5ErrorManager {
   public:
      /// HDF5ErrorManager class' default constructor
      HDF5ErrorManager();

      /// HDF5ErrorManager class' default constructor
      ~HDF5ErrorManager();

      /// Checks the existance of the main HDF5 file (data.h5)
      static void hdf5_file_checker();

      /// Checks the hdf5 dataset layout
      static void hdf5_dataset_layout_checker(std::string dataset_name, H5D_layout_t layout);
};

}

#endif // HDF5_ERROR_MANAGER_H

