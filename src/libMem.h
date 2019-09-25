#ifndef LIBMEM_H
#define LIBMEM_H

namespace mcpdft {

class LibMem {
   public:
      /// Query memory information from linux
      void query_system_memory(const struct sysinfo *info);
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
#endif // LIBMEM_H
