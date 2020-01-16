#ifndef IWRITE_H
#define IWRITE_H

#include <armadillo>

namespace mcpdft {

class IWrite {
   public:
      /// Calls write_opdm() and write_tpdm() fxns
      virtual void write_rdms(const arma::mat &D1a,
                              const arma::mat &D1b,
			      const arma::mat &D2ab) = 0;

      /// Write 1RDM into memory/disk
      virtual void write_opdm(const arma::mat &D1a,
		              const arma::mat &D1b) = 0;

      /// Write 2RDM into memory/disk 
      virtual void write_tpdm(const arma::mat &D2ab) = 0;
};

}

#endif // IWRITE_H
