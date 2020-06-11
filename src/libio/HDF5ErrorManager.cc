#include <string>
#include <fstream>
#include "HDF5ErrorManager.h"
#include "TOC.h"

namespace mcpdft {

   HDF5ErrorManager::HDF5ErrorManager() {};

   HDF5ErrorManager::~HDF5ErrorManager() {};

   void HDF5ErrorManager::hdf5_file_checker() {
      std::ifstream file(H5FILE);
      if( !file.good() ) {
         throw "\n  Warning: No accessible HDF5 file by the name \"data.h5\" exists!\n";
         file.close();
      }
   }

   void HDF5ErrorManager::hdf5_dataset_layout_checker(std::string dataset_name, H5D_layout_t layout) {
      std::printf ("  Storage layout for %s is   : ", dataset_name.c_str());
      switch (layout) {
          case H5D_COMPACT:
              printf ("H5D_COMPACT\n");
              break;
          case H5D_CONTIGUOUS:
              printf ("H5D_CONTIGUOUS\n");
              break;
          case H5D_CHUNKED:
              printf ("H5D_CHUNKED\n");
              break;
          case H5D_VIRTUAL:
              printf ("H5D_VIRTUAL\n");
              break;
          case H5D_LAYOUT_ERROR:
          case H5D_NLAYOUTS:
              printf ("H5D_LAYOUT_ERROR\n");
      }
   }
}
