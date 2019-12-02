#ifndef IREAD_H
#define IREAD_H

#include <armadillo>

namespace mcpdft {

class IRead {
   public:
      /// Calls read_opdm() and read_tpdm() fxns
      virtual void read_rdms(arma::mat &D1a,
                             arma::mat &D1b,
			     arma::mat &D2ab) = 0;

      /// Read 1RDM into memory/disk
      virtual void read_opdm(arma::mat &D1a,
		             arma::mat &D1b) = 0;

      /// Read 1RDM into memory/disk
      virtual void read_tpdm(arma::mat &D2ab) = 0;
};

}

#endif // IREAD_H
