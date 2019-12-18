#ifndef HDF5COMPACTFACTORY_H
#define HDF5COMPACTFACTORY_H

#include "IOFactory.h"

namespace mcpdft {

class HDF5CompactFactory: public IOFactory {
   public:
      IRead* create_IRead() override {
         return new HDF5ReadCompact;
      }

      IWrite* create_IWrite() override {
         return new HDF5WriteCompact;
      }
};

}

#endif   // HDF5COMPACTFACTORY_H

