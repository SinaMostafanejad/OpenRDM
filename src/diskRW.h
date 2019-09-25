#ifndef DISKRW_H
#define DISKRW_H

namespace mcpdft {

    class DiskRW: public IRDMInterface {
        public:
	   /// The default Constructor
           DiskRW(); 		

           /// Destructor
	   ~DiskRW();

           /// Calculates the required amount of memory needed for dealing with RDMs
           virtual void calculate_memory(arma::mat &D1, arma::mat &D2ab);

           /// Calls read_opdm() and read_tpdm() fxns
           virtual void read_rdms();

           /// Read 1RDM into disk
           virtual void read_opdm();

           /// Read 1RDM into disk
           virtual void read_tpdm();

           /// Calls write_opdm() and write_tpdm() fxns
           virtual void write_rdms();

           /// Write 1RDM into disk
           virtual void write_opdm();

           /// Write 2RDM into disk
           virtual void write_tpdm();
	private:
	   /// A memory buffer
           long int memory_;
};


}

#endif // DISKRW_H
