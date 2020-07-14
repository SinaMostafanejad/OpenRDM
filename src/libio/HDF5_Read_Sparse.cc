#include "hdf5.h"
#include "TOC.h"
#include "HDF5Utility.h"

namespace mcpdft {

   void HDF5Utility::read_sparse_coo_opdm(arma::vec &val,
					  arma::Col<int> &row_idx,
					  arma::Col<int> &col_idx,
		                          const bool is_active,
					  const std::string &rdm_type) {

      /* file indentifiers and handles */
      hid_t file_id;
//      hid_t nnz_dst_id;
      hid_t val_dst_id, row_idx_dst_id, col_idx_dst_id;
      herr_t status;

      /* open the existing HDF5 file in the read-only mode */
      file_id = H5Fopen(H5FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* open the existing HDF5 dataset in the read-only mode */
      if (rdm_type == "D1A") {
         if (is_active) { /* is_ao == false && is_active == true */
//            nnz_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1A_MO_NNZ, H5P_DEFAULT);
            val_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1A_MO_VAL, H5P_DEFAULT);
            row_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1A_MO_ROW_IDX, H5P_DEFAULT);
            col_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1A_MO_COL_IDX, H5P_DEFAULT);
         }else{ /* is_ao == false && is_active == false */
//            nnz_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1A_MO_NNZ, H5P_DEFAULT);
            val_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1A_MO_VAL, H5P_DEFAULT);
            row_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1A_MO_ROW_IDX, H5P_DEFAULT);
            col_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1A_MO_COL_IDX, H5P_DEFAULT);
         }
      } else if (rdm_type == "D1B") {
         if (is_active) { /* is_ao == false && is_active == true */
//            nnz_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1B_MO_NNZ, H5P_DEFAULT);
            val_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1B_MO_VAL, H5P_DEFAULT);
            row_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1B_MO_ROW_IDX, H5P_DEFAULT);
            col_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_FULL_D1B_MO_COL_IDX, H5P_DEFAULT);
         }else{ /* is_ao == false && is_active == false */
//            nnz_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1B_MO_NNZ, H5P_DEFAULT);
            val_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1B_MO_VAL, H5P_DEFAULT);
            row_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1B_MO_ROW_IDX, H5P_DEFAULT);
            col_idx_dst_id = H5Dopen2(file_id, H5D_SP_SYM_ACT_D1B_MO_COL_IDX, H5P_DEFAULT);
         }
      } else {
         throw "\n Error: rdm_type should match \"D1A\" or \"D1B\".\n";
      }
      /* read the dataSet */
//      status = H5Dread(nnz_dst_id, H5T_NATIVE_INT,
//                       H5S_ALL, H5S_ALL, H5P_DEFAULT, &nnz);

      status = H5Dread(val_dst_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, val.memptr());

      status = H5Dread(row_idx_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, row_idx.memptr());

      status = H5Dread(col_idx_dst_id, H5T_NATIVE_INT,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT, col_idx.memptr());

      /* close datasets */
//      status = H5Dclose(nnz_dst_id);
      status = H5Dclose(val_dst_id);
      status = H5Dclose(row_idx_dst_id);
      status = H5Dclose(col_idx_dst_id);

      /* close the file */
      status = H5Fclose(file_id);
   }

}
