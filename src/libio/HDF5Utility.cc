#include "hdf5.h"
#include "HDF5ErrorManager.h"
#include "HDF5Utility.h"
#include "TOC.h"

namespace mcpdft {

   void HDF5Utility::read_nbfs() {
      HDF5ErrorManager::hdf5_file_checker();
      int *nao{nullptr}, *nmo{nullptr}; 

      /* file indentifiers and handles */
      hid_t file_id;
      hid_t nao_dst_id, nmo_dst_id;
      herr_t status;

      /* Retrieve the data set name from TOC */
      std::string nao_dst_name = H5D_N_AO;
      std::string nmo_dst_name = H5D_N_MO;

      /* Open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      nao_dst_id = H5Dopen2(file_id, H5D_N_AO, H5P_DEFAULT);
      nmo_dst_id = H5Dopen2(file_id, H5D_N_MO, H5P_DEFAULT);

      /* Read the DataSet */
      status = H5Dread(nao_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, nao);

      status = H5Dread(nmo_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, nmo);

      /* Close datasets */
      status = H5Dclose(nao_dst_id);
      status = H5Dclose(nmo_dst_id);

      /* Close the file */
      status = H5Fclose(file_id);
   }
   
   void HDF5Utility::read_grids() {
      HDF5ErrorManager::hdf5_file_checker();
   }

}
