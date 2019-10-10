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
   
       size_t npts = mc->get_npts();
       arma::vec W(mc->get_w());

       double exc = 0.0;
       for (int p = 0; p < npts; p++) {
           double rhoa = rho_a(p);
           double rhob = rho_b(p);
           double rhoa_43 = pow( rhoa, 4.0/3.0); 
           double rhob_43 = pow( rhob, 4.0/3.0); 
           double Xa = sqrt(sigma_aa(p)) / rhoa_43;
           double Xb = sqrt(sigma_bb(p)) / rhob_43;
           
           double Sa = (Xa * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
           double Sb = (Xb * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
           
           double Fsa = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sa,2.0)) / KAPPA ), -1.0 );
           double Fsb = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sb,2.0)) / KAPPA ), -1.0 );
          
           auto E = [](double rhos, double Fss) -> double{
                    double temp = -0.75 * pow(3.0, 1.0/3.0) * pow(M_PI, 2.0/3.0) * pow(rhos,4.0/3.0) * Fss / M_PI;
                    return temp;
           };

           double EX_GGAa = 0.5 * E(2.0*rhoa,Fsa);
           double EX_GGAb = 0.5 * E(2.0*rhob,Fsb);
          
           exc += ( EX_GGAa + EX_GGAb ) * W(p); 
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

       double exc = 0.0;
       for (int p = 0; p < npts; p++) {
           double rhoa = rho_a(p);
           double rhob = rho_b(p);
   
           double sigmaaa = sigma_aa(p);
           double sigmaab = sigma_ab(p);
           double sigmabb = sigma_bb(p);
   
           double rho = rhoa + rhob;
           double sigma = sigmaaa + sigmabb + 2.0 * sigmaab;
   
           if ( rho > tol ) {
              if ( rhoa < tol ){
                 rho = rhob;
                 sigmabb = std::max(0.0,sigmabb);
                 sigma = sigmabb;
                 double t2 = 1.0 / rhob;
                 double t3 = pow(t2 ,1.0/3.0);
                 double t6 = pow(t2 ,1.0/6.0);
                 double t9 = sqrt(t2);
                 double t11 = t3 * t3;
                 double t17 = log(1.0 + 0.3216395899738507e2 / (0.1112037486309468e2 * t6 + 0.3844746237447211e1 * t3 + 0.1644733775567609e1
                            * t9 + 0.2405871291288192 * t11));
                 double t18 = (1.0 + 0.1274696188700087 * t3) * t17;
                 double t20 = pow(rhob ,2.0);
                 double t21 = pow(rhob ,1.0/3.0);
                 double t23 = 1.0 / t21 / t20;
                 double t26 = exp(0.2000000587336264e1 * t18);
                 double t27 = t26 - 1.0;
                 double t31 = 0.2162211495206379 / t27 * sigmabb * t23;
                 double t33 = pow(t27 ,2.0);
                 double t35 = pow(sigmabb ,2.0);
                 double t37 = pow(t20 ,2.0);
                 double t38 = pow(t21 ,2.0);
                 double t49 = log(1.0 + 0.2162211495206379 * sigmabb * t23 * (1.0 + t31) /(1.0 + t31 + 0.4675158550002605e-1 / t33 * t35 / t38 / t37));
                 double zk = rhob * (-0.310907e-1 * t18 + 0.1554534543482745e-1 * t49);
   
                 exc += rho * zk * W(p);
              }else if ( rhob < tol ){
                       rho = rhoa;
                       sigmaaa = std::max(0.0,sigmaaa);
                       sigma = sigmaaa;
                       double t2 = 1.0 / rhoa;
                       double t3 = pow(t2 ,1.0/3.0);
                       double t6 = pow(t2 ,1.0/6.0);
                       double t9 = sqrt(t2);
                       double t11 = t3 * t3;
                       double t17 = log(1.0 + 0.3216395899738507e2 / (0.1112037486309468e2 * t6 + 0.3844746237447211e1 * t3 + 0.1644733775567609e1
                                  * t9 + 0.2405871291288192 * t11));
                       double t18 = (1.0 + 0.1274696188700087 * t3) * t17;
                       double t20 = pow(rhoa ,2.0);
                       double t21 = pow(rhoa ,1.0/3.0);
                       double t23 = 1.0 / t21 / t20;
                       double t26 = exp(0.2000000587336264e1 * t18);
                       double t27 = t26 - 1.0;
                       double t31 = 0.2162211495206379 / t27 * sigmaaa * t23;
                       double t33 = pow(t27 ,2.0);
                       double t35 = pow(sigmaaa ,2.0);
                       double t37 = pow(t20 ,2.0);
                       double t38 = pow(t21 ,2.0);
                       double t49 = log(1.0 + 0.2162211495206379 * sigmaaa * t23 * (1.0 + t31) /(1.0 + t31 + 0.4675158550002605e-1 / t33 * t35 / t38 / t37));
                       double zk = rhoa * (-0.310907e-1 * t18 + 0.1554534543482745e-1 * t49);
   
                       exc += rho * zk * W(p);
              }else{
                   double t4 = 1/rho;
                   double t5 = pow( t4, 1.0/3.0);
                   double t7 = 1.0 + 0.1325688999052018 * t5;
                   double t8 = pow(t4, 1.0/6.0);
                   double t11 = sqrt(t4);
                   double t13 = pow(t5 ,2.0);
                   double t15 = 0.598255043577108e1 * t8 + 0.2225569421150687e1 * t5 + 0.8004286349993634 * t11 + 0.1897004325747559 * t13;
                   double t18 = 1.0 + 0.1608197949869254e2 / t15;
                   double t19 = log(t18);
                   double t21 = 0.621814e-1 * t7 * t19;
                   double t23 = 1.0 + 0.6901399211255825e-1 * t5;
                   double t28 = 0.8157414703487641e1 * t8 + 0.2247591863577616e1 * t5 + 0.4300972471276643 * t11 + 0.1911512595127338 * t13;
                   double t31 = 1.0 + 0.2960874997779344e2 / t28;
                   double t32 = log(t31);
                   double t33 = t23 * t32;
                   double t35 = rhoa - 1.0 * rhob;
                   double t36 = t35 * t4;  // zeta
                   double t37 = 1.0 + t36;
                   double t38 = pow(t37 ,1.0/3.0);
                   double t41 = 1.0 - t36;
                   double t42 = pow(t41 ,1.0/3.0);
                   double t44 = t38 * t37 + t42 * t41 - 2.0;
                   double t45 = pow(t35 ,2.0);
                   double t46 = pow(t45 ,2.0);
                   double t47 = pow(rho ,2.0);
                   double t48 = pow(t47 ,2.0);
                   double t49 = 1.0 / t48;
                   double t50 = t46 * t49;
                   double t52 = 1.0 - t50;
                   double t55 = 0.37995525e-1 * t33 * t44 * t52;
                   double t57 = 1.0 + 0.1274696188700087 * t5;
                   double t62 = 0.1112037486309468e2 * t8 + 0.3844746237447211e1 * t5 + 0.1644733775567609e1 * t11 + 0.2405871291288192 * t13;
                   double t65 = 1.0 + 0.3216395899738507e2 / t62;
                   double t66 = log(t65);
                   double t69 = -0.310907e-1 * t57 * t66 + t21;
                   double t70 = t69 * t44;
                   double t72 = 0.1923661050931536e1 * t70 * t50;
                   double t73 = pow(t38 ,2.0);
                   double t75 = pow(t42 ,2.0);
                   double t77 = 0.5 * t73 + 0.5 * t75;
                   double t78 = pow(t77 ,2.0);
                   double t79 = t78 * t77;
                   double t80 = 1.0 / t78;
                   double t81 = sigma * t80;
                   double t82 = pow(rho ,1.0/3.0);
                   double t84 = 1.0 / t82 / t47;
                   double t85 = -t21 + t55 + t72;
                   double t86 = 1.0 / t79;
                   double t89 = exp(-0.3216396844291482e2 * t85 * t86);
                   double t90 = t89 - 1.0;
                   double t91 = 1.0 / t90;
                   double t92 = t91 * sigma;
                   double t93 = t80 * t84;
                   double t95 = 0.1362107888567592 * t92 * t93;
                   double t96 = 1.0 + t95;
                   double t98 = pow(t90 ,2.0);
                   double t99 = 1.0 / t98;
                   double t100 = pow(sigma ,2.0);
                   double t101 = t99 * t100;
                   double t102 = pow(t78 ,2.0);
                   double t103 = 1.0 / t102;
                   double t104 = pow(t82 ,2.0);
                   double t106 = 1.0 / t104 / t48;
                   double t107 = t103 * t106;
                   double t110 = 1.0 + t95 + 0.1855337900098064e-1 * t101 * t107;
                   double t111 = 1.0 / t110;
                   double t115 = 1.0 + 0.1362107888567592 * t81 * t84 * t96 * t111;
                   double t116 = log(t115);
                   double t118 = 0.310906908696549e-1 * t79 * t116;
                   double zk = -t21 + t55 + t72 + t118;
   
                   exc += rho * zk * W(p);
              }
              }else{
                   double zk = 0.0;
                   exc += rho * zk * W(p);
                   }
       }
       return exc;
   }

}
