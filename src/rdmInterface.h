#ifndef RDMINTERFACE_H
#define RDMINTERFACE_H

#include <armadillo>

namespace mcpdft {

class IRDMInterface {
   public:
      /// Virtual destructor in case we want to delete an IRDMInterface pointer
      virtual ~IRDMInterface() {}

      /// Calculates the required amount of memory needed for dealing with RDMs
      virtual void calculate_memory(arma::mat &D1, arma::mat &D2ab) = 0;
      
      /// Calls read_opdm() and read_tpdm() fxns
      virtual void read_rdms() = 0;

      /// Read 1RDM into memory/disk
      virtual void read_opdm() = 0;

      /// Read 1RDM into memory/disk
      virtual void read_tpdm() = 0;

      /// Calls write_opdm() and write_tpdm() fxns
      virtual void write_rdms() = 0;

      /// Write 1RDM into memory/disk
      virtual void write_opdm() = 0;

      /// Write 2RDM into memory/disk 
      virtual void write_tpdm() = 0;
};

}

#endif // RDMINTERFACE_H
