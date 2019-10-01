#include <string>
#include <stdio.h>
#include <fstream>
#include <limits>
#include <sys/sysinfo.h>
#include "libMem.h"

namespace mcpdft {
   void LibMem::query_system_memory(const struct sysinfo *info) {
      mem_unit_   = info->mem_unit;
      totalram_   = info->totalram;
      freeram_    = info->freeram;
      sharedram_  = info->sharedram;
      bufferram_  = info->bufferram;
      totalswap_  = info->totalswap;
      freeswap_   = info->freeswap;

      unsigned long totalram_MB;
      unsigned long freeram_MB;
      unsigned long sharedram_MB;
      unsigned long bufferram_MB;
      unsigned long totalswap_MB;
      unsigned long freeswap_MB;

      totalram_MB   = (mem_unit_ * totalram_  ) /1024 /1024;
      freeram_MB    = (mem_unit_ * freeram_   ) /1024 /1024;
      sharedram_MB  = (mem_unit_ * sharedram_ ) /1024 /1024;
      bufferram_MB  = (mem_unit_ * bufferram_ ) /1024 /1024;
      totalswap_MB  = (mem_unit_ * totalswap_ ) /1024 /1024;
      freeswap_MB   = (mem_unit_ * freeswap_  ) /1024 /1024;

      printf("\n===================================================\n");
      printf("   Total usable main memory size = %-8lu (MB)\n",totalram_MB);
      printf("   Available memory size         = %-8lu (MB)\n",freeram_MB);
      printf("   Amount of shared memory       = %-8lu (MB)\n",sharedram_MB);
      printf("   Buffer memory size            = %-8lu (MB)\n",bufferram_MB);
      printf("   Total swap space size         = %-8lu (MB)\n",totalswap_MB);
      printf("   Available swap sapce          = %-8lu (MB)\n",freeswap_MB);
      printf("===================================================\n");
   }
   unsigned long LibMem::get_mem_total() {
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
