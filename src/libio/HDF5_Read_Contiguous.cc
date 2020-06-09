#include <armadillo>
#include <fstream>
#include "mcpdft.h"
#include "HDF5_Read_Contiguous.h"
#include "hdf5.h"
#include <assert.h>

#include "TOC.h"

#define RANK 2

namespace mcpdft {

   void HDF5ReadContiguous::read_rdms(arma::mat &D1a,
                                      arma::mat &D1b,
			              arma::mat &D2ab) {
      read_opdm(D1a, D1b);
      read_tpdm(D2ab);
   }

   void HDF5ReadContiguous::read_opdm(arma::mat &D1a,
		                      arma::mat &D1b) {
      std::ifstream file(OPDM_H5FILE);
      if( !file.good() )
	throw "\n  Warning: No accessible HDF5 file by the name \"opdm.h5\" exists!\n";
        file.close();

      assert( D1a.n_cols == D1a.n_rows );
      assert( D1b.n_cols == D1b.n_rows );
      size_t dim = D1a.n_cols;

      /* file indentifiers and handles */
      hid_t file_id;
      hid_t dcpl_D1a, dcpl_D1b;
      hid_t D1a_dst_id, D1b_dst_id;
      herr_t status;

      /* Open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(OPDM_H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      D1a_dst_id = H5Dopen2(file_id, D1A_DATASET, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      D1b_dst_id = H5Dopen2(file_id, D1B_DATASET, H5P_DEFAULT);

      /*
       * Retrieve the dataset creation property list for opdms
       * and print their storage layout.
       */
      dcpl_D1a = H5Dget_create_plist (D1a_dst_id);
      dcpl_D1b = H5Dget_create_plist (D1b_dst_id);
      H5D_layout_t layout[] = { H5Pget_layout (dcpl_D1a), H5Pget_layout (dcpl_D1b) };
      std::string  dataset_names[] = { D1A_DATASET, D1B_DATASET };

      printf ("-------------------------------------------------------------\n");
      printf ("          HDF5 OPDM & TPDM FILES LAYOUT SUMMARY              \n");
      printf ("-------------------------------------------------------------\n");
      for (int i = 0; i < 2; i++) {
         printf ("  Storage layout for %s is   : ", dataset_names[i].c_str());
         switch (layout[i]) {
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

      /* Read the D1a_DataSet */
      status = H5Dread(D1a_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, D1a.memptr());

      /* Read the D1a_DataSet */
      status = H5Dread(D1b_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, D1b.memptr());

      /* Close property lists */
      status = H5Pclose(dcpl_D1a);
      status = H5Pclose(dcpl_D1b);

      /* Close datasets */
      status = H5Dclose(D1a_dst_id);
      status = H5Dclose(D1b_dst_id);

      /* Close the file */
      status = H5Fclose(file_id); 
   }

   void HDF5ReadContiguous::read_tpdm(arma::mat &D2ab) {
      std::ifstream file(TPDM_H5FILE);
      if( !file.good() )
	throw "\n  Warning: No accessible HDF5 file by the name \"tpdm.h5\" exists!\n";
        file.close();

      assert( D2ab.n_cols == D2ab.n_rows );
      size_t dim = D2ab.n_cols;

      /* file indentifiers and handles */
      hid_t file_id;
      hid_t dcpl_D2ab;
      hid_t D2ab_dst_id;
      herr_t status;

      /* Open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(TPDM_H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* Open the existing HDF5 dataset in the read-only mode */
      D2ab_dst_id = H5Dopen2(file_id, D2AB_DATASET, H5P_DEFAULT);

      /*
       * Retrieve the dataset creation property list for opdms
       * and print their storage layout
       */
      dcpl_D2ab = H5Dget_create_plist (D2ab_dst_id);
      H5D_layout_t layout = H5Pget_layout (dcpl_D2ab);

      printf ("  Storage layout for %s is : ", D2AB_DATASET);
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
      printf ("-------------------------------------------------------------\n\n" );

      /* Read the D2ab_DataSet */
      status = H5Dread(D2ab_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, D2ab.memptr());

      /* Close property list */
      status = H5Pclose(dcpl_D2ab);
      
      /* Close dataset */
      status = H5Dclose(D2ab_dst_id);

      /* Close the file */
      status = H5Fclose(file_id); 
   }
}
