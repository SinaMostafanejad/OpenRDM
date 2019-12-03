#ifndef HDF5_READ_H
#define HDF5_READ_H

#include <armadillo>
#include "IRead.h"

namespace mcpdft {

class HDF5_Read: public IRead {
    public:
       /// Calls read_opdm() and read_tpdm() fxns
       void read_rdms(arma::mat &D1a,
                      arma::mat &D1b,
                      arma::mat &D2ab);

       /// Read 1RDM from disk
       void read_opdm(arma::mat &D1a,
    	              arma::mat &D1b);

       /// Read 2RDM from disk
       void read_tpdm(arma::mat &D2ab);
};

}

#endif // HDF5_READ_H
