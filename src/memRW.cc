#include <armadillo>
#include <sys/sysinfo.h>
#include <stdio.h>
#include "memRW.h"

namespace mcpdft {
    MemRW::MemRW()  {}
    MemRW::~MemRW() {}

    void MemRW::calculate_memory(arma::mat &D1, arma::mat &D2ab) {
       struct sysinfo info;
       sysinfo(&info);

       mem_unit_ = info.mem_unit;
       totalram_ = info.totalram;
       freeram_  = info.freeram;

       unsigned long totalram_MB;
       unsigned long freeram_MB;
       totalram_MB = (mem_unit_ * totalram_) /1024 /1024;
       freeram_MB  = (mem_unit_ * freeram_ ) /1024 /1024;

       printf("\n===================================================\n");
       printf("   Total usable main memory size:  %-8lu (MB)\n",totalram_MB);
       printf("   Available memory size:          %-8lu (MB)\n",freeram_MB);
       printf("===================================================\n");
    }

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
