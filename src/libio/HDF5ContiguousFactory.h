#ifndef HDF5CONTIGUOUSFACTORY_H
#define HDF5CONTIGUOUSFACTORY_H

#include "IOFactory.h"

namespace mcpdft {

class HDF5ContiguousFactory: public IOFactory {
   public:
      IRead* create_IRead() override {
         return new HDF5ReadContiguous;
      }

      IWrite* create_IWrite() override {
         return new HDF5WriteContiguous;
      }
};

}

#endif   // HDF5CONTIGUOUSFACTORY_H

