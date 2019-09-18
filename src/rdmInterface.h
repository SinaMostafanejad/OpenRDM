#ifndef RDMINTERFACE_H
#define RDMINTERFACE_H

#include <armadillo>
#include "mcpdft.h"

namespace mcpdft {

class IRDMInterface {

   public:
       
      /// Reads the 1- and 2- RDMs into memory
      void virtual read_into_memory() = 0;

      /// Reads the 1- and 2-RDMs into disk batches
      void virtual read_into_disk() = 0;

      /// Calculate the required amount of memory needed for dealing with RDMs
      void virtual calculate_memory(arma::mat &D1, arma::mat &D2ab) = 0;

      /// Write 1RDM into disk
      void virtual write_opdm_into_disk() = 0;

      /// Write 2RDM into dist 
      void virtual write_tpdm_into_disk() = 0;
 
};

}

#endif // RDMINTERFACE_H
