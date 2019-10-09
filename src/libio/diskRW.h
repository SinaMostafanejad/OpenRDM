#ifndef DISKRW_H
#define DISKRW_H

#include <armadillo>
#include "rdmInterface.h"
#include "mcpdft.h"

namespace mcpdft {

    class DiskRW: public IRDMInterface {
        public:
	   /// The default Constructor
           DiskRW(); 		

           /// Destructor
	   ~DiskRW();

//           /// Calculates the required amount of memory needed for dealing with RDMs
//           virtual void calculate_memory(arma::mat &D1, arma::mat &D2ab);

           /// Calls read_opdm() and read_tpdm() fxns
           void read_rdms();

           /// Read 1RDM into disk
           void read_opdm(arma::mat &D1a,
		          arma::mat &D1b);

           /// Read 1RDM into disk
           void read_tpdm();

           /// Calls write_opdm() and write_tpdm() fxns
           void write_rdms();

           /// Write 1RDM into disk
           void write_opdm(const arma::mat &D1a,
			   const arma::mat &D1b);

           /// Write 2RDM into disk
           void write_tpdm();
	private:
	   /// A memory buffer
           long int memory_;
};


}

#endif // DISKRW_H
