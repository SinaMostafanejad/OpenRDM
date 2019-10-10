#include <armadillo>
#include "mcpdft.h"
#include "functional.h"

namespace mcpdft {
   /*================================================================*/
   /*                      Exchange Functionals                      */
   /*================================================================*/
   Functional::Functional()  {}
   Functional::~Functional() {}

   double Functional::EX_LSDA(const MCPDFT *mc,
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

   double Functional::EX_PBE(const MCPDFT *mc,
		             const arma::vec &rho_a,
		             const arma::vec &rho_b,
		             const arma::vec &sigma_aa,
		             const arma::vec &sigma_bb) {

       const double delta = 0.06672455060314922;
       const double MU = (1.0/3.0) * delta * M_PI * M_PI;
       const double KAPPA = 0.804;
       double tol = 1.0e-20;
   
       auto kF = [](double RHO) -> double {
                 double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
                 return dum;
       };
   
       auto eX = [=](double RHO) -> double {
                 double temp = -(3.0 * kF(RHO)) / (4.0 * M_PI);
                 return temp;
       };
   
       auto FX = [=](double SIGMA) -> double {
                 double temp = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(SIGMA,2.0)) / KAPPA ), -1.0 );
                 return temp;
       };
   
   
       auto S = [=](double RHO, double SIGMA) -> double {
                double temp = sqrt(SIGMA) / (2.0 * kF(RHO) * RHO);
                return temp;
       };
   
       size_t npts = mc->get_npts();
       arma::vec W(mc->get_w());

       double exc = 0.0;
       for (int p = 0; p < npts; p++) {
           double rhoa = rho_a(p);
           double rhob = rho_b(p);
           double rho = rhoa + rhob;
           double sigmaaa = sigma_aa(p);
           double sigmabb = sigma_bb(p);
           double sigma = 0.0;
           if ( rho > tol ) {
              if ( rhoa < tol ){
                 rho = rhob;
                 sigmabb = std::max(0.0,sigmabb);
                 sigma = sigmabb;
   
                 double zk = eX(2.0 * rho) * FX(S(2.0 * rho, 4.0 * sigma));
                 exc += rho * zk * W(p);
              }else if ( rhob < tol ){
                       rho = rhoa;
                       sigmaaa = std::max(0.0,sigmaaa);
                       sigma = sigmaaa;
                       double zk = eX(2.0 * rho) * FX(S(2.0 * rho, 4.0 * sigma));
                       exc += rho * zk * W(p);
              }else {
                    double zka = rhoa * eX(2.0 * rhoa) * FX(S(2.0 * rhoa, 4.0 * sigmaaa));
                    double zkb = rhob * eX(2.0 * rhob) * FX(S(2.0 * rhob, 4.0 * sigmabb));
                    double zk = zka + zkb;
                    exc += zk * W(p);
              }
           }else{
                   exc += 0.0;
                }
       }
       return exc;
   }
   /*================================================================*/
   /*                    Correlation Functionals                     */
   /*================================================================*/
    double Functional::EC_VWN3(const MCPDFT *mc,
                               const arma::vec &rho_a,
                               const arma::vec &rho_b) {
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

   double Functional::EC_PBE(const MCPDFT *mc,
		             const arma::vec &rho_a,
		             const arma::vec &rho_b,
                             const arma::vec &sigma_aa,
                             const arma::vec &sigma_ab,
			     const arma::vec &sigma_bb) {
      double tol = 1.0e-20;
      size_t npts = mc->get_npts();
      arma::vec W(mc->get_w());
      const double pa = 1.0;
      const double Aa = 0.0168869;
      const double a1a = 0.11125;
      const double b1a = 10.357;
      const double b2a = 3.6231;
      const double b3a = 0.88026;
      const double b4a = 0.49671;
      const double pe = 1.0;
      const double c0p = 0.0310907;
      const double a1p = 0.21370;
      const double b1p = 7.5957;
      const double b2p = 3.5876;
      const double b3p = 1.6382;
      const double b4p = 0.49294;
      const double c0f = 0.01554535;
      const double a1f = 0.20548;
      const double b1f = 14.1189;
      const double b2f = 6.1977;
      const double b3f = 3.3662;
      const double b4f = 0.62517;
      const double d2Fz = 1.7099209341613656173;
      const double BETA = 0.06672455060314922;
      const double GAMMA = 0.0310906908696549;

      auto Fi = [=](double ZETA) -> double {
                double dumm = 0.5 * (pow((1.0 + ZETA) ,2.0/3.0 ) + pow((1.0 - ZETA) ,2.0/3.0));
                return dumm;
      };

      auto kF = [](double RHO) -> double {
                double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
                return dum;
      };

      auto ks = [=](double RHO) -> double {
                double temp = sqrt(4.0 * kF(RHO) / M_PI);
                return temp;
      };

      auto Fz = [](double ZETA) -> double {
                double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
                return dum;
      };

      auto t = [=](double RHO, double SIGMA, double ZETA) -> double {
               double temp = sqrt(SIGMA) / (2.0 * ks(RHO) * Fi(ZETA) * RHO);
               return temp;
      };

      auto G = [](double r, double T, double a1, double b1, double b2, double b3, double b4, double p) -> double {
               double dum = -2.0 * T * (1.0 + a1 * r) * log(1.0 + 0.5 * pow(T * (b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0)) ,-1.0));
               return dum;

      };

      auto Ac = [=](double r) -> double {
                double temp = -G(r,Aa,a1a,b1a,b2a,b3a,b4a,pa);
                return temp;
      };

      auto EcP = [=](double r) -> double {
                 double dum = G(r,c0p,a1p,b1p,b2p,b3p,b4p,pe);
                 return dum;
      };

      auto EcF = [=](double r) -> double {
                 double dumm = G(r,c0f,a1f,b1f,b2f,b3f,b4f,pe);
                 return dumm;
      };

      auto Ec = [=](double r, double ZETA) -> double {
                double dum = EcP(r) + ( Ac(r) * Fz(ZETA) * (1.0 - pow(ZETA ,4.0)) ) / d2Fz + ( EcF(r) - EcP(r) ) * Fz(ZETA) * pow(ZETA ,4.0);
                return dum;
      };

      auto A = [=](double r, double ZETA) -> double {
               double dum = (BETA/GAMMA) * pow( exp(-Ec(r,ZETA) / (pow(Fi(ZETA),3.0) * GAMMA)) - 1.0, -1.0);
               return dum;
      };

      auto H = [=](double RHO, double SIGMA, double r, double ZETA) -> double {
               double temp = pow(Fi(ZETA),3.0) * GAMMA * log(1.0 + (BETA/GAMMA) * pow(t(RHO,SIGMA,ZETA) ,2.0) * (1.0 + A(r,ZETA) * pow(t(RHO,SIGMA,ZETA),2.0))
                           / (1.0 + A(r,ZETA) * pow(t(RHO,SIGMA,ZETA),2.0) + pow(A(r,ZETA),2.0) * pow(t(RHO,SIGMA,ZETA),4.0)));
               return temp;
      };

      double exc = 0.0;
      for (int p = 0; p < npts; p++) {
          double rhoa = rho_a(p);
          double rhob = rho_b(p);
          double rho = rhoa + rhob;
          double zeta = (rhoa - rhob) / rho;
          double rs =  pow( 3.0 / ( 4.0 * M_PI * rho) , 1.0/3.0 );
          double sigmaaa = sigma_aa(p);
          double sigmaab = sigma_ab(p);
          double sigmabb = sigma_bb(p);
          double sigma = sigmaaa + sigmabb + 2.0 * sigmaab;
          if ( rho > tol ) {
             if ( rhoa < tol ){
                rho = rhob;
                sigmabb = std::max(0.0,sigmabb);
                sigma = sigmabb;
                zeta = 1.0;
             }else if ( rhob < tol ){
                      rho = rhoa;
                      sigmaaa = std::max(0.0,sigmaaa);
                      sigma = sigmaaa;
                      zeta = 1.0;
             }else/* if (!(rhoa < tol) && !(rhob < tol) ) */{
                      sigmaaa = std::max(0.0,sigmaaa);
                      sigmabb = std::max(0.0,sigmabb);
                      sigma = sigmaaa + sigmabb + 2.0 * sigmaab;
             }
             double zk = H(rho,sigma,rs,zeta) + Ec(rs,zeta);
             exc += rho * zk * W(p);
          }else{
                  double zk = 0.0;
                  exc += rho * zk * W(p);
               }
      }
      return exc;
   }
}
