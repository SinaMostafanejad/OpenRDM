#ifndef RDMINTERFACE_H
#define RDMINTERFACE_H

#include <armadillo>

namespace mcpdft {

class IRDMInterface {
   public:
      /// Virtual destructor in case we want to delete an IRDMInterface pointer
      virtual ~IRDMInterface() {}

      /// Calls read_opdm() and read_tpdm() fxns
      virtual void read_rdms(arma::mat &D1a,
                             arma::mat &D1b,
			     arma::mat &D2ab) = 0;

      /// Read 1RDM into memory/disk
      virtual void read_opdm(arma::mat &D1a,
		             arma::mat &D1b) = 0;

      /// Read 1RDM into memory/disk
      virtual void read_tpdm(arma::mat &D2ab) = 0;

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

#endif // RDMINTERFACE_H
