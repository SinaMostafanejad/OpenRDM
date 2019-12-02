#include <armadillo>
#include <fstream>
#include "mcpdft.h"
#include "HDF5_Write.h"
#include "hdf5.h"
#include <assert.h>

#define OPDM_H5FILE  "opdm.h5"
#define TPDM_H5FILE  "tpdm.h5"
#define D1A_DATASET  "/D1a/D1a_DataSet"
#define D1B_DATASET  "/D1b/D1b_DataSet"
#define D2AB_DATASET "/D2ab/D2ab_DataSet"
#define RANK 2

namespace mcpdft {
   void HDF5_Write::write_rdms(const arma::mat &D1a, 
                               const arma::mat &D1b,
			       const arma::mat &D2ab) {
      write_opdm(D1a, D1b);
      write_tpdm(D2ab);
   }

   void HDF5_Write::write_opdm(const arma::mat &D1a,
		               const arma::mat &D1b) {
      assert( D1a.n_cols == D1a.n_rows );
      assert( D1b.n_cols == D1b.n_rows );
      size_t dim = D1a.n_cols;

      /* file indentifiers and handles */
      hid_t file_id;
      hid_t D1a_grp_id, D1b_grp_id;
      hid_t D1a_dst_id, D1b_dst_id;
      hid_t dataspace_id;
      hid_t filespace_id;//, memspace_id;
      hid_t prop_id;
      /* dataset dimensions at creation time */
      double D1_buff[dim][dim];
      hsize_t dims[2] = {dim, dim};
      // hsize_t      maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
      herr_t       status;
      // hsize_t      chunk_dims[2] = {2, 5};

      /* Create a new file. If the file exists, its contents
       * will be overwritten */
      file_id = H5Fcreate(OPDM_H5FILE, H5F_ACC_TRUNC,
	                  H5P_DEFAULT, H5P_DEFAULT);

      /* Create "D1a" and "D1b" groups within the root group */
      D1a_grp_id = H5Gcreate2(file_id, "/D1a",
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      D1b_grp_id = H5Gcreate2(file_id, "/D1b",
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /* Create the data space for the datasets */
      dataspace_id = H5Screate_simple(RANK, dims, NULL);

      /* Create a dataset in group "/D1a" */
      D1a_dst_id = H5Dcreate2(file_id, D1A_DATASET,
        	              H5T_IEEE_F64LE, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /* Create a dataset in group "/D1b" */
      D1b_dst_id = H5Dcreate2(file_id, D1B_DATASET,
        	              H5T_IEEE_F64LE, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      for(int i=0; i < dim; i++)
         for(int j=0; j < dim; j++)
            D1_buff[i][j] = D1a(i,j);

      /* Write the D1a_DataSet */
      status = H5Dwrite(D1a_dst_id, H5T_NATIVE_DOUBLE,
        	        H5S_ALL, H5S_ALL, H5P_DEFAULT, D1_buff);

      for(int i=0; i < dim; i++)
         for(int j=0; j < dim; j++)
            D1_buff[i][j] = D1b(i,j);

      /* Write the D1b_DataSet */
      status = H5Dwrite(D1b_dst_id, H5T_NATIVE_DOUBLE,
        	        H5S_ALL, H5S_ALL, H5P_DEFAULT, D1_buff);

      // /* Modify dataset creation properties, i.e. enable chunking  */
      // prop_id = H5Pcreate (H5P_DATASET_CREATE);
      // status = H5Pset_chunk (prop_id, RANK, chunk_dims);

      /* Close dataspace. */
      status = H5Sclose(dataspace_id);

      /* Close datasets. */
      status = H5Dclose(D1a_dst_id);
      status = H5Dclose(D1b_dst_id);

      /* Close groups. */
      status = H5Gclose(D1a_grp_id);
      status = H5Gclose(D1b_grp_id);

      /* Close the file. */
      status = H5Fclose(file_id); 
   }

   void HDF5_Write::write_tpdm(const arma::mat &D2ab) {
      assert( D2ab.n_cols == D2ab.n_rows );
      size_t dim = D2ab.n_cols;

      /* file indentifiers and handles */
      hid_t file_id;
      hid_t D2ab_grp_id;
      hid_t D2ab_dst_id;
      hid_t dataspace_id;
      hid_t filespace_id;//, memspace_id;
      hid_t prop_id;
      /* dataset dimensions at creation time */
      double D2_buff[dim][dim];
      hsize_t dims[2] = {dim, dim};
      // hsize_t      maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
      herr_t       status;
      // hsize_t      chunk_dims[2] = {2, 5};

      /* Create a new file. If the file exists, its contents
       * will be overwritten */
      file_id = H5Fcreate(TPDM_H5FILE, H5F_ACC_TRUNC,
	                  H5P_DEFAULT, H5P_DEFAULT);

      /* Create "D2ab" group within the root group */
      D2ab_grp_id = H5Gcreate2(file_id, "/D2ab",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /* Create the data space for the dataset */
      dataspace_id = H5Screate_simple(RANK, dims, NULL);

      /* Create a dataset in group "/D2ab" */
      D2ab_dst_id = H5Dcreate2(file_id, D2AB_DATASET,
        	              H5T_IEEE_F64LE, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      for(int ij=0; ij < dim; ij++)
         for(int kl=0; kl < dim; kl++)
             D2_buff[ij][kl] = D2ab(ij,kl);

      /* Write the D1a_DataSet */
      status = H5Dwrite(D2ab_dst_id, H5T_NATIVE_DOUBLE,
        	        H5S_ALL, H5S_ALL, H5P_DEFAULT, D2_buff);

      // /* Modify dataset creation properties, i.e. enable chunking  */
      // prop_id = H5Pcreate (H5P_DATASET_CREATE);
      // status = H5Pset_chunk (prop_id, RANK, chunk_dims);

      /* Close dataspace. */
      status = H5Sclose(dataspace_id);

      /* Close dataset. */
      status = H5Dclose(D2ab_dst_id);

      /* Close group. */
      status = H5Gclose(D2ab_grp_id);

      /* Close the file. */
      status = H5Fclose(file_id); 
   }
}
