#include "HDF5Client.h"
#include "HDF5Utility.h"

namespace mcpdft {

   HDF5Client::HDF5Client() {};

   HDF5Client::~HDF5Client() {};

   void HDF5Client::factory_client(H5D_layout_t layout,
                                   factory_mode mode,
				   arma::mat &D1a,
                                   arma::mat &D1b,
                                   arma::mat &D2ab) {
      IOFactory* iof;
      IRead*  ird;
      IWrite* iwt;
  
      switch (layout) {
        case H5D_COMPACT:
           iof = new HDF5CompactFactory;
           break;
        case H5D_CONTIGUOUS:
           iof = new HDF5ContiguousFactory;
           break;
        case H5D_CHUNKED:
           iof = new HDF5ChunkedFactory;
           break;
        case H5D_VIRTUAL:
           //printf ("H5D_VIRTUAL\n");
           //break;
        case H5D_LAYOUT_ERROR:
        case H5D_NLAYOUTS:
           throw "Warning: H5D_LAYOUT_ERROR!\n";
      }

      switch(mode) {
         case READ:
	    {	 
               ird = iof->create_IRead();
               ird->read_rdms(D1a,D1b,D2ab);
	       break;
	    }
	 case WRITE:
	    {
               iwt = iof->create_IWrite();
               // iwt->write_rdms(d1a,d1b,d2ab);
	       break;
	    }
      }

      delete iof;
   }

#if 0
   arma::mat HDF5Client::get_D1a() const { return D1a_;} 
   arma::mat HDF5Client::get_D1b() const { return D1b_;} 
   arma::mat HDF5Client::get_D2ab() const { return D2ab_;} 

   void HDF5Client::set_D1a(const arma::mat &D1a)  { D1a_  = D1a;} 
   void HDF5Client::set_D1b(const arma::mat &D1b)  { D1b_  = D1b;} 
   void HDF5Client::set_D2ab(const arma::mat &D2ab) { D2ab_ = D2ab;} 
#endif

}
