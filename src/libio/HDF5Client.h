#ifndef HDF5CLIENT_H
#define HDF5CLIENT_H

#include <armadillo>
#include <string>
//#include "IOFactory.h"
#include "HDF5_Read_Contiguous.h"
#include "HDF5_Write_Contiguous.h"
#include "HDF5_Read_Compact.h"
#include "HDF5_Write_Compact.h"
#include "HDF5_Read_Chunked.h"
#include "HDF5_Write_Chunked.h"
#include "HDF5ContiguousFactory.h"
#include "HDF5CompactFactory.h"
#include "HDF5ChunkedFactory.h"

namespace mcpdft {

class HDF5Client {
   public:
      /// constructor
      HDF5Client();

      /// destructor
      ~HDF5Client();

      enum factory_mode {
         READ ,
	 WRITE,
      };

      /// HDF5 factory client read
      void factory_client(H5D_layout_t layout,
		          factory_mode mode,
                          const arma::mat &D1a,
	                  const arma::mat &D1b,
	                  const arma::mat &D2ab);

//      factory_client_write(H5D_layout_t layout)

      /* accessors */
      void set_factory_layout(const H5D_layout_t layout);

      H5D_layout_t get_factory_layout(H5D_layout_t layout) const;

      arma::mat get_D1a() const;
      arma::mat get_D1b() const;
      arma::mat get_D2ab() const;

      void set_D1a(const arma::mat &D1a);
      void set_D1b(const arma::mat &D1b);
      void set_D2ab(const arma::mat &D2ab);

   private:
//       /// IOFactory factory object
//       IOFactory* io_factory_;
//       
//       /// IRead factory object
//       IRead* iread_factory_;
// 
//       /// IWrite factory object
//       IWrite* iwrite_factory;

      /// HDF5 factory (layout) type
      H5D_layout_t layout_;

      /// factory mode (R and W or RW modes)
      std::string factory_mode_;

      /// alpha-spin 1-electron reduced-density matrix
      arma::mat D1a_;

      /// alpha-spin 1-electron reduced-density matrix
      arma::mat D1b_;

      /// alpha-beta spin-block of 2-electron reduced-density matrix
      arma::mat D2ab_;

};

}

#endif // HDF5CLIENT_H
