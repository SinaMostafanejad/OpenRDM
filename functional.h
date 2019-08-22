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
      double EX_LSDA(arma::vec &rho_a,
                     arma::vec &rho_b,
                     MCPDFT* mc);

      /// VWN RPA expression-3 correlation functional
      double EC_VWN3_RPA_III(arma::vec &rho_a,
                             arma::vec &rho_b,
                             MCPDFT* mc);

};

}

#endif // FUNCTIONAL_H 
