#include <armadillo>
#include "mcpdft.h"
#include "diskRW.h"
#include "hdf5.h"
#include <assert.h>

#define OPDM_H5FILE "opdm.h5"
#define TPDM_H5FILE "tpdm.h5"

namespace mcpdft {
   DiskRW::DiskRW()  {}
   DiskRW::~DiskRW() {}

   void DiskRW::read_rdms() {
   }

   void DiskRW::read_opdm(const arma::mat &D1a,
		          const arma::mat &D1b) {
      assert( D1a.n_cols == D1a.n_rows );
      assert( D1b.n_cols == D1b.n_rows );
      size_t dim = D1a.n_cols;

      /* file indentifiers and handles */
      hid_t file;
      hid_t dataspace, dataset;
      hid_t filespace, memspace;
      hid_t prop;
      /* dataset dimensions at creation time */
      hsize_t dims[2] = {dim, dim};
      hsize_t      maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
      herr_t       status;
      hsize_t      chunk_dims[2] = {2, 5};
      int          data[3][3];
   }

   void DiskRW::read_tpdm() {}

   void DiskRW::write_rdms() {}

   void DiskRW::write_opdm() {}

   void DiskRW::write_tpdm() {}
}
