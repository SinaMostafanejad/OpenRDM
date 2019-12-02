#ifndef IOFACTORY_H
#define IOFACTORY_H

#include "IRead.h"
#include "IWrite.h"

namespace mcpdft {

class IOFactory {
   public:
      virtual IRead*  create_IRead()  = 0;
      virtual IWrite* create_IWrite() = 0;
};

}

#endif   // IOFACTORY_H
