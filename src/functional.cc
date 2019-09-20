#include <armadillo>
#include "mcpdft.h"
#include "functional.h"


namespace mcpdft {

   Functional::Functional()  {}
   Functional::~Functional() {}

   double Functional::EX_LSDA(const MCPDFT* mc,
                              const arma::vec &rho_a,
                              const arma::vec &rho_b) {
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

    double Functional::EC_VWN3(const MCPDFT* mc,
                               const arma::vec& rho_a,
                               const arma::vec& rho_b) {
       double tol = 1.0e-20;

       const double ecp1 = 0.03109070000;
       const double ecp2 = -0.409286;
       const double ecp3 = 13.0720;
       const double ecp4 = 42.7198;
       const double ecf1 = 0.01554535000;
       const double ecf2 = -0.743294;
       const double ecf3 = 20.1231;
       const double ecf4 = 101.578;
       const double d2Fz = 1.7099209341613656173;

       size_t npts = mc->get_npts();
       arma::vec W(mc->get_w());

       auto x = [](double RHO) -> double {
                double rs = pow( 3.0 / ( 4.0 * M_PI * RHO ) , 1.0/3.0 );
                double dum = sqrt(rs);
                return dum;
       };
       auto Fz = [](double ZETA) -> double {
                 double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
                 return dum;
       };
       auto X = [](double i, double c, double d) -> double {
                double temp = pow(i,2.0) + c * i + d;
                return temp;
       };
       auto Q = [](double c, double d) -> double {
                double temp1 = sqrt( 4 * d - pow(c,2.0) );
                return temp1;
       };
       auto q = [=](double RHO, double A, double p, double c, double d) -> double {
                double dum1 = A * ( log( pow(x(RHO),2.0) / X(x(RHO),c,d) ) + 2.0 * c * atan( Q(c,d)/(2.0*x(RHO) + c) ) * pow(Q(c,d),-1.0)
                            - c * p * ( log( pow(x(RHO)-p,2.0) / X(x(RHO),c,d) ) + 2.0 * (c + 2.0 * p) * atan( Q(c,d)/(2.0*x(RHO) + c) )
                            * pow(Q(c,d),-1.0) ) * pow(X(p,c,d),-1.0) );
                return dum1;
       };
       auto EcP = [=](double RHO) -> double {
                  double dumm = q(RHO,ecp1,ecp2,ecp3,ecp4);
                  return dumm;
       };
       auto EcF = [=](double RHO) -> double {
                  double dum = q(RHO,ecf1,ecf2,ecf3,ecf4);
                  return dum;
       };

       double exc = 0.0;
       for (int p = 0; p < npts; p++) {
           double rhoa = rho_a(p);
           double rhob = rho_b(p);
           double rho = rhoa + rhob;
           double zeta = 0.0;
           if ( rho > tol ) {
              if ( rhoa < tol ){
                 rho = rhob;
                 zeta = 1.0;
              }else if ( rhob < tol ){
                       rho = rhoa;
                       zeta = 1.0;
              }else {/* if (!(rhoa < tol) && !(rhob < tol) ) */
                    zeta = (rhoa - rhob) / rho;
              }
              double zk = EcP(rho) + Fz(zeta) * (EcF(rho) - EcP(rho));
              exc += rho * zk * W(p);
           }else{
                   double zk = 0.0;
                   exc += 0.0;
                }
       }
       return exc;
   }

}
