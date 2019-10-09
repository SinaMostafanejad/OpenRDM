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
      double EX_LSDA(const MCPDFT *mc,
                     const arma::vec &rho_a,
                     const arma::vec &rho_b);

      /// Perdew-Burke-Ernzerhof (PBE) exchange functional
      double EX_PBE(const MCPDFT *mc,
                    const arma::vec &rho_a,
                    const arma::vec &rho_b,
                    const arma::vec &sigma_aa,
                    const arma::vec &sigma_bb);

      /*================================================================*/
      /*                    Correlation Functionals                     */
      /*================================================================*/
      /// VWN RPA expression-3 correlation functional
      double EC_VWN3(const MCPDFT *mc,
                     const arma::vec &rho_a,
                     const arma::vec &rho_b);

      /// Perdew-Burke-Ernzerhof (PBE) correlation functional
      double EC_PBE(const MCPDFT *mc,
                    const arma::vec &rho_a,
                    const arma::vec &rho_b,
                    const arma::vec &sigma_aa,
                    const arma::vec &sigma_ab,
                    const arma::vec &sigma_bb);
};

}

#endif // FUNCTIONAL_H 
