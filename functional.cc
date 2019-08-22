#include <armadillo>
#include "mcpdft.h"
#include "functional.h"


namespace mcpdft {

   Functional::Functional()  {};
   Functional::~Functional() {};

   double Functional::EX_LSDA(arma::vec &rho_a,
                              arma::vec &rho_b,
                              MCPDFT* mc) {

       const double alpha = (2.0/3.0);      // Slater value: a constant
       const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);

       size_t npts = mc->get_npts();
       arma::vec W(mc->get_w());

       double exc = 0.0;
       for (size_t p = 0; p < npts; p++) {
           double exa = pow(2.0,1.0/3.0) * Cx * pow( rho_a(p), 4.0/3.0) ;
           double exb = pow(2.0,1.0/3.0) * Cx * pow( rho_b(p), 4.0/3.0) ;
           double ex_LSDA = exa + exb;
           exc += - ex_LSDA * W(p);
       }
       return exc;
   }

}
