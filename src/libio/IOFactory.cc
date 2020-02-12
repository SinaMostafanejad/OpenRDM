#include "IOFactory.h"
#include "HDF5Factory.h"

namespace mcpdft {
   IOFactory*  IOFactory::createFactory(HDF5Layout layout) {
      switch(layout) {
         case CONTIGUOUS:
            return new HDF5Factory;
//	 case COMPACT:
//            return new HDF5CompactFactory();
//	 case CHUNKED:
//            return new HDF5ChunkedFactory();
//	 case VIRTUAL:
//            return new HDF5VirtualFactory();
	 default:
            printf ("H5D_LAYOUT_ERROR\n");
	    break;
      }
}
