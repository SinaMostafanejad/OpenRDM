#ifndef HDF5CHUNKEDFACTORY_H
#define HDF5CHUNKEDFACTORY_H

#include "IOFactory.h"

namespace mcpdft {

class HDF5ChunkedFactory: public IOFactory {
   public:
      IRead* create_IRead() override {
         return new HDF5ReadChunked;
      }

      IWrite* create_IWrite() override {
         return new HDF5WriteChunked;
      }
};

}

#endif   // HDF5CHUNKEDFACTORY_H

