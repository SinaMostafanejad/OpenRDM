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
};

}
#endif // MEMRW_H
