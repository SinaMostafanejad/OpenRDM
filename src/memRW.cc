#include <armadillo>
#include <stdio.h>
#include "memRW.h"

namespace mcpdft {
    MemRW::MemRW()  {}
    MemRW::~MemRW() {}

    void MemRW::read_rdms() {}

    void MemRW::read_opdm() {}

    void MemRW::read_tpdm() {}

    void MemRW::write_rdms() {}

    void MemRW::write_opdm() {}

    void MemRW::write_tpdm() {}

    unsigned long MemRW::get_mem_total() {
       std::string token;
       std::ifstream file("/proc/meminfo");
       while(file >> token) {
           if(token == "MemTotal:") {
               unsigned long mem;
               if(file >> mem) {
                   return mem;
                   file.close();
               } else {
                   return 0;
               }
           }
           // ignore rest of the line
           file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
           file.close();
       }
       return 0; // nothing found
    }
}
