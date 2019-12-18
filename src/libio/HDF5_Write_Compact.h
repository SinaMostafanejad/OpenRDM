#ifndef HDF5_WRITE_COMPACT_H
#define HDF5_WRITE_COMPACT_H

#include <armadillo>
#include "IWrite.h"

namespace mcpdft {

    class HDF5WriteCompact: public IWrite {
        public:
           /// Calls write_opdm() and write_tpdm() fxns
           void write_rdms(const arma::mat &D1a,
                           const arma::mat &D1b,
			   const arma::mat &D2ab);

           /// Write 1RDM into disk
           void write_opdm(const arma::mat &D1a,
			   const arma::mat &D1b);

           /// Write 2RDM into disk
           void write_tpdm(const arma::mat &D2ab);
};


}

#endif // HDF5_WRITE_COMPACT_H
