#include "hdf5.h"
#include "HDF5ErrorManager.h"
#include "HDF5Utility.h"
#include "TOC.h"

namespace mcpdft {

   HDF5Utility::HDF5Utility() {
      HDF5ErrorManager::hdf5_file_checker();
   }

   HDF5Utility::~HDF5Utility() {}

   void HDF5Utility::read_nbfs(size_t &nao,
		               size_t &nmo,
			       size_t &npts) {
      /* file indentifiers and handles */
      hid_t file_id;
      hid_t nao_dst_id, nmo_dst_id, npts_dst_id;
      herr_t status;

      /* Open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      nao_dst_id  = H5Dopen2(file_id, H5D_N_AO, H5P_DEFAULT);
      nmo_dst_id  = H5Dopen2(file_id, H5D_N_MO, H5P_DEFAULT);
      npts_dst_id = H5Dopen2(file_id, H5D_N_PTS, H5P_DEFAULT);

      /* Read the DataSet */
      status = H5Dread(nao_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nao);

      status = H5Dread(nmo_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nmo);

      status = H5Dread(npts_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &npts);

      /* Close datasets */
      status = H5Dclose(nao_dst_id);
      status = H5Dclose(nmo_dst_id);
      status = H5Dclose(npts_dst_id);

      /* Close the file */
      status = H5Fclose(file_id);
   }

   void HDF5Utility::read_active_space_vars(size_t &nactele,
		                            size_t &nactorb,
			                    size_t &ncore, 
			                    size_t &nfrz, 
			                    size_t &nocc, 
			                    size_t &nvir) {
      /* file indentifiers and handles */
      hid_t file_id;
      hid_t nactele_dst_id, nactorb_dst_id,
      ncore_dst_id, nfrz_dst_id, nocc_dst_id, nvir_dst_id;
      herr_t status;

      /* Open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      nactele_dst_id = H5Dopen2(file_id, H5D_N_CAS_ELE, H5P_DEFAULT);
      nactorb_dst_id = H5Dopen2(file_id, H5D_N_CAS_ORB, H5P_DEFAULT);
      ncore_dst_id   = H5Dopen2(file_id, H5D_N_COR, H5P_DEFAULT);
      nfrz_dst_id    = H5Dopen2(file_id, H5D_N_FRZ, H5P_DEFAULT);
      nocc_dst_id    = H5Dopen2(file_id, H5D_N_OCC_ACT, H5P_DEFAULT);
      nvir_dst_id    = H5Dopen2(file_id, H5D_N_VIR_ACT, H5P_DEFAULT);

      /* Read the DataSet */
      status = H5Dread(nactele_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nactele);

      status = H5Dread(nactorb_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nactorb);

      status = H5Dread(ncore_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &ncore);

      status = H5Dread(nfrz_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nfrz);

      status = H5Dread(nocc_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nocc);

      status = H5Dread(nvir_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nvir);

      /* Close datasets */
      status = H5Dclose(nactele_dst_id);
      status = H5Dclose(nactorb_dst_id);
      status = H5Dclose(ncore_dst_id);
      status = H5Dclose(nfrz_dst_id);
      status = H5Dclose(nocc_dst_id);
      status = H5Dclose(nvir_dst_id);

      /* Close the file */
      status = H5Fclose(file_id);
   }

   void HDF5Utility::read_AO2MO_Cmat(arma::mat &Cmat) {
      /* file indentifiers and handles */
      hid_t file_id;
      hid_t cmat_dst_id;
      herr_t status;

      /* open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* open the existing HDF5 dataset in the read-only mode */
      cmat_dst_id = H5Dopen2(file_id, H5D_C, H5P_DEFAULT);

      /* read the dataSet */
      status = H5Dread(cmat_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, Cmat.memptr());

      /* close datasets */
      status = H5Dclose(cmat_dst_id);

      /* close the file */
      status = H5Fclose(file_id);
   }

   void HDF5Utility::read_grids(arma::vec &W,
		                arma::vec &X,
				arma::vec &Y,
				arma::vec &Z) {
      /* file indentifiers and handles */
      hid_t file_id;
      hid_t w_dst_id, x_dst_id, y_dst_id, z_dst_id;
      herr_t status;

      /* open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* open the existing HDF5 dataset in the read-only mode */
      w_dst_id = H5Dopen2(file_id, H5D_W, H5P_DEFAULT);
      x_dst_id = H5Dopen2(file_id, H5D_X, H5P_DEFAULT);
      y_dst_id = H5Dopen2(file_id, H5D_Y, H5P_DEFAULT);
      z_dst_id = H5Dopen2(file_id, H5D_Z, H5P_DEFAULT);

      /* read the dataSet */
      status = H5Dread(w_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, W.memptr());

      status = H5Dread(x_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, X.memptr());

      status = H5Dread(y_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, Y.memptr());

      status = H5Dread(z_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, Z.memptr());

      /* check to see if the number of points match */
//       if( W.n_elem != X.n_elem != Y.n_elem != Z.n_elem )
//          throw "\n Warning: The sizes of grid vectors w, x, y and z do not match.\n";

      /* close datasets */
      status = H5Dclose(w_dst_id);
      status = H5Dclose(x_dst_id);
      status = H5Dclose(y_dst_id);
      status = H5Dclose(z_dst_id);

      /* close the file */
      status = H5Fclose(file_id);
   }

}
