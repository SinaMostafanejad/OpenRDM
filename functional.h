#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include <armadillo>
#include "mcpdft.h"

namespace mcpdft {

class Functional {

   public:
 
      /// constructor
      Functional();

      /// destructor
      ~Functional();

      /// Slater exchange functional 
      double EX_LSDA(MCPDFT* mc,
                     arma::vec &rho_a,
                     arma::vec &rho_b);

      /// VWN RPA expression-3 correlation functional
      double EC_VWN3(MCPDFT* mc,
                     arma::vec &rho_a,
                     arma::vec &rho_b);

};

}

#endif // FUNCTIONAL_H 
