#ifndef MEMRW_H
#define MEMRW_H

#include "rdmInterface.h"

namespace mcpdft {

class MemRW : public IRDMInterface {
    public:
       /// The default Constructor
       MemRW(); 		

       /// Destructor
       ~MemRW();

       /// Calculates the required amount of memory needed for dealing with RDMs
       virtual void calculate_memory(arma::mat &D1, arma::mat &D2ab);

       /// Calls read_opdm() and read_tpdm() fxns
       virtual void read_rdms();

       /// Read 1RDM into memory
       virtual void read_opdm();

       /// Read 1RDM into memory
       virtual void read_tpdm();

       /// Calls write_opdm() and write_tpdm() fxns
       virtual void write_rdms();

       /// Write 1RDM into memory
       virtual void write_opdm();

       /// Write 2RDM into memory
       virtual void write_tpdm();

       /// Get total memory using C++ function
       unsigned long get_mem_total();

    private:
       /// Uptime run since boot in seconds
       long uptime_;

       /// 1, 5 and 15 minutes' load averages
       unsigned long loads_[3];

       /// Total usable main memory size
       unsigned long totalram_;

       /// Available memory size
       unsigned long freeram_;

       /// Amount of shared memory
       unsigned long sharedram_;

       /// Buffer memory
       unsigned long bufferram_;

       /// Total swap space size
       unsigned long totalswap_;

       /// Available swap space
       unsigned long freeswap_;

       /// Number of ongoing processes
       unsigned short num_procs_;

       /// Total high memory size
       unsigned long totalhigh_;

       /// Available high memory size
       unsigned long freehigh_;

       /// Memory unit size in bytes
       unsigned int mem_unit_;

       /// Padding for libc5 
       char _f[20-2*sizeof(long)-sizeof(int)];
};

}
#endif // MEMRW_H
