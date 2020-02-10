#ifndef HDF5CLIENT_H
#define HDF5CLIENT_H

#include <armadillo>
#include <string>
#include "HDF5_Read_Contiguous.h"
#include "HDF5_Write_Contiguous.h"
#include "HDF5_Read_Compact.h"
#include "HDF5_Write_Compact.h"
#include "HDF5_Read_Chunked.h"
#include "HDF5_Write_Chunked.h"
#include "HDF5ContiguousFactory.h"
#include "HDF5CompactFactory.h"
#include "HDF5ChunkedFactory.h"
#include "hdf5.h"

namespace mcpdft {

class HDF5Client {
   public:
      /// constructor
      HDF5Client();

      /// destructor
      ~HDF5Client();

      /// factory modes
      enum factory_mode {
         READ ,
	 WRITE,
      };

      /// HDF5 factory client
      void factory_client(H5D_layout_t layout,
		          factory_mode mode,
                          const arma::mat &D1a,
	                  const arma::mat &D1b,
	                  const arma::mat &D2ab);

      /* accessors */

      /// setting out factory layout
      void set_factory_layout(const H5D_layout_t layout);
 
      /// getting factory layout
      H5D_layout_t get_factory_layout() const;

      /// get alpha-spin 1-electron reduced-density matrix
      arma::mat get_D1a() const;
      /// get beta-spin 1-electron reduced-density matrix
      arma::mat get_D1b() const;
      /// get alpha-beta spin-block of the 2-electron reduced-density matrix
      arma::mat get_D2ab() const;

      /// set alpha-spin 1-electron reduced-density matrix
      void set_D1a(const arma::mat &D1a);
      /// set beta-spin 1-electron reduced-density matrix
      void set_D1b(const arma::mat &D1b);
      /// set alpha-beta spin-block of the 2-electron reduced-density matrix
      void set_D2ab(const arma::mat &D2ab);

   private:
      /// HDF5 factory (layout) type
      H5D_layout_t layout_;

      /// factory mode (READ and WRITE)
      std::string factory_mode_;

      /// alpha-spin 1-electron reduced-density matrix
      arma::mat D1a_;

      /// alpha-spin 1-electron reduced-density matrix
      arma::mat D1b_;

      /// alpha-beta spin-block of the 2-electron reduced-density matrix
      arma::mat D2ab_;

};

}

#endif // HDF5CLIENT_H
