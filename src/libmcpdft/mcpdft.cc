#include <armadillo>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include "mcpdft.h"
#include "openrdmConfig.h"

#ifdef WITH_OPENMP // _OPENMP
   #include <omp.h>
#endif

namespace mcpdft {

   MCPDFT::MCPDFT(std::string test_case)  { common_init(test_case); }
   MCPDFT::MCPDFT() {}
   MCPDFT::~MCPDFT() {}

   void MCPDFT::common_init(std::string test_case) {
       print_banner();
       read_grids_from_file(test_case);
       read_orbitals_from_file(test_case);
       read_energies_from_file(test_case);
       if(test_case == "h2_tpbe_sto3g") {
	  is_gga_ = true;
	  read_gradients_from_file(test_case);
       }else{
	  is_gga_ = false;
       }

       read_opdm_from_file(test_case);
       // read_cmat_from_file();
   }

   void MCPDFT::build_rho() {
      build_density_functions();
      if ( is_gga_ ) {
         build_density_gradients();
      }
   }

   void MCPDFT::build_density_functions() {
      int nbfs = get_nbfs();
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());
      arma::vec W(get_w());
      arma::vec rhoa(npts, arma::fill::zeros);
      arma::vec rhob(npts, arma::fill::zeros);
      arma::vec rho(npts, arma::fill::zeros);
      double dum_a = 0.0;
      double dum_b = 0.0;
      double dum_tot = 0.0;
      size_t chunk_size = 0;
      int p{0}, mu{0}, nu{0};
      #ifdef WITH_OPENMP
         int nthrds{0};
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
         printf("                   *** Warning ***\n");
         printf("   Calculating the density (gradients) using OpenMP");
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
         nthrds = omp_get_max_threads();
         nthrds /= 2;
         omp_set_num_threads(nthrds);
      #endif

      #pragma omp parallel default(shared) \
                           private(p, mu, nu)
      {  
         #pragma omp for schedule(static) \
                         reduction(+:dum_a, dum_b, dum_tot) \
	                 nowait
            for(p = 0; p < npts; p++) {
               double tempa = 0.0;
               double tempb = 0.0;
               #pragma omp parallel for schedule(static) \
                                        reduction(+:tempa,tempb) \
	                                num_threads(2) \
                                        collapse(2)
                  for(mu = 0; mu < nbfs; mu++) {
                     for(int nu = 0; nu < nbfs; nu++) {
                        tempa += D1a(mu, nu) * phi(p, mu) * phi(p, nu);
                        tempb += D1b(mu, nu) * phi(p, mu) * phi(p, nu);
                     }
                  }
                  rhoa(p) = tempa;
                  rhob(p) = tempb;
                  rho(p) = rhoa(p) + rhob(p);

                  dum_a += rhoa(p) * W(p);
                  dum_b += rhob(p) * W(p);
                  dum_tot += ( rhoa(p) + rhob(p) ) * W(p) ;
            } /* end of omp parallel for loop */
      } /* end of omp parallel region */
      set_rhoa(rhoa);
      set_rhob(rhob);
      set_rho(rho);

      printf("\n");
      printf("  Integrated total density = %20.12lf\n",dum_tot);
      printf("  Integrated alpha density = %20.12lf\n",dum_a);
      printf("  Integrated beta density  = %20.12lf\n",dum_b);
      printf("\n");
   }

   void MCPDFT::build_density_gradients() {
      int nbfs = get_nbfs();
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());
      int p{0}, nu{0}, sigma{0};
      arma::mat phi_x(get_phi_x());
      arma::mat phi_y(get_phi_y());
      arma::mat phi_z(get_phi_z());
      arma::vec rho_a_x(npts, arma::fill::zeros);
      arma::vec rho_b_x(npts, arma::fill::zeros);
      arma::vec rho_a_y(npts, arma::fill::zeros);
      arma::vec rho_b_y(npts, arma::fill::zeros);
      arma::vec rho_a_z(npts, arma::fill::zeros);
      arma::vec rho_b_z(npts, arma::fill::zeros);
      arma::vec sigma_aa(npts, arma::fill::zeros);
      arma::vec sigma_ab(npts, arma::fill::zeros);
      arma::vec sigma_bb(npts, arma::fill::zeros);
      #pragma omp parallel default(shared) \
                           private(p, nu, sigma)
      {
         #pragma omp for schedule(static)
             for (int p = 0; p < npts; p++) {
                 double duma_x = 0.0;
                 double dumb_x = 0.0;
                 double duma_y = 0.0;
                 double dumb_y = 0.0;
                 double duma_z = 0.0;
                 double dumb_z = 0.0;
                 #pragma omp parallel for schedule(static) \
                   	                  reduction(+:duma_x, duma_y, duma_z,\
                   		                      dumb_x, dumb_y, dumb_z)\
                                          shared(p, nbfs, \
    				                 sigma_aa, sigma_bb, sigma_ab, \
    				                 rho_a_x, rho_a_y, rho_a_z, \
    				                 rho_b_x, rho_b_y, rho_b_z) \
    			                  num_threads(2) \
                                          collapse(2)
                    for (int sigma = 0; sigma < nbfs; sigma++) {
                        for (int nu = 0; nu < nbfs; nu++) {
                            duma_x += ( phi_x(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_x(p, nu) ) * D1a(sigma, nu);
                            dumb_x += ( phi_x(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_x(p, nu) ) * D1b(sigma, nu);
                            duma_y += ( phi_y(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_y(p, nu) ) * D1a(sigma, nu);
                            dumb_y += ( phi_y(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_y(p, nu) ) * D1b(sigma, nu);
                            duma_z += ( phi_z(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_z(p, nu) ) * D1a(sigma, nu);
                            dumb_z += ( phi_z(p, sigma) * phi(p, nu) + phi(p, sigma) * phi_z(p, nu) ) * D1b(sigma, nu);
                        }
                    }
                    rho_a_x(p) = duma_x;
                    rho_b_x(p) = dumb_x;
                    rho_a_y(p) = duma_y;
                    rho_b_y(p) = dumb_y;
                    rho_a_z(p) = duma_z;
                    rho_b_z(p) = dumb_z;
                    sigma_aa(p) = ( rho_a_x(p) * rho_a_x(p) ) +  ( rho_a_y(p) * rho_a_y(p) ) + ( rho_a_z(p) * rho_a_z(p) );
                    sigma_bb(p) = ( rho_b_x(p) * rho_b_x(p) ) +  ( rho_b_y(p) * rho_b_y(p) ) + ( rho_b_z(p) * rho_b_z(p) );
                    sigma_ab(p) = ( rho_a_x(p) * rho_b_x(p) ) +  ( rho_a_y(p) * rho_b_y(p) ) + ( rho_a_z(p) * rho_b_z(p) );
             }
      }
      set_rhoa_x(rho_a_x);
      set_rhob_x(rho_b_x);
      set_rhoa_y(rho_a_y);
      set_rhob_y(rho_b_y);
      set_rhoa_z(rho_a_z);
      set_rhob_z(rho_b_z);
      set_sigma_aa(sigma_aa);
      set_sigma_ab(sigma_ab);
      set_sigma_bb(sigma_bb);
   }

   void MCPDFT::build_pi(const arma::mat &D2ab) {
      build_ontop_pair_density(D2ab);
      if ( is_gga_ ) {
         build_ontop_pair_density_gradients(D2ab);
      }
   }

   void MCPDFT::build_ontop_pair_density(const arma::mat &D2ab) {
      int nbfs = get_nbfs();
      size_t npts = get_npts();
      arma::vec temp(npts);
      arma::mat phi(get_phi());
      int p{0},
	  mu{0}, nu{0},
	  lambda{0}, sigma{0};
      for (int p = 0; p < npts; p++) {
          double dum = 0.0;
          // pi(r,r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
          #pragma omp parallel for default(shared) \
          	                       private(mu,nu,lambda,sigma) \
          	                       reduction(+:dum) \
          	                       collapse(4)
             for (int mu = 0; mu < nbfs; mu++) {
                 for (int nu = 0; nu < nbfs; nu++) {
                     for (int lambda = 0; lambda < nbfs; lambda++) {
                         for (int sigma = 0; sigma < nbfs; sigma++) {
                             dum += phi(p, mu) * phi(p, nu) * phi(p, lambda) * phi(p, sigma) * D2ab(nu*nbfs+mu, sigma*nbfs+lambda);
                         }
                     }
                 }
             } /* end of omp parallel for loop */
             temp(p) = dum;
      }
      set_pi(temp);
   }

   void MCPDFT::build_ontop_pair_density_gradients(const arma::mat &D2ab) {
      int nbfs = get_nbfs();
      size_t npts = get_npts();
      int p{0},
	  mu{0}, nu{0},
	  lambda{0}, sigma{0};
      arma::mat phi(get_phi());
      arma::mat phi_x(get_phi_x());
      arma::mat phi_y(get_phi_y());
      arma::mat phi_z(get_phi_z());
      arma::vec pi_x(npts, arma::fill::zeros);
      arma::vec pi_y(npts, arma::fill::zeros);
      arma::vec pi_z(npts, arma::fill::zeros);
      #pragma omp parallel default(shared) \
                           private(p, mu, nu, lambda, sigma)
      {
         #pragma omp for schedule(static)
            for (int p = 0; p < npts; p++) {
                double dum_x = 0.0;
                double dum_y = 0.0;
                double dum_z = 0.0;
                // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
                #pragma omp parallel for schedule(static) \
                  	                 reduction(+:dum_x, dum_y, dum_z) \
    		                         num_threads(2) \
                                         collapse(4)
                for (int mu = 0; mu < nbfs; mu++) {
                    for (int nu = 0; nu < nbfs; nu++) {
                        for (int lambda = 0; lambda < nbfs; lambda++) {
                            for (int sigma = 0; sigma < nbfs; sigma++) {
                                dum_x += ( phi_x(p, mu) * phi(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi_x(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi_x(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi(p, sigma) * phi_x(p, nu) ) * D2ab(nu*nbfs+mu, sigma*nbfs+lambda);

                                dum_y += ( phi_y(p, mu) * phi(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi_y(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi_y(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi(p, sigma) * phi_y(p, nu) ) * D2ab(nu*nbfs+mu, sigma*nbfs+lambda);

                                dum_z += ( phi_z(p, mu) * phi(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi_z(p, lambda) * phi(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi_z(p, sigma) * phi(p, nu) +
                                           phi(p, mu) * phi(p, lambda) * phi(p, sigma) * phi_z(p, nu) ) * D2ab(nu*nbfs+mu, sigma*nbfs+lambda);
                            }
                        }
                    }
                }
                pi_x(p) = dum_x;
                pi_y(p) = dum_y;
                pi_z(p) = dum_z;
            }
      }
      set_pi_x(pi_x);
      set_pi_z(pi_y);
      set_pi_z(pi_z);
   }

   void MCPDFT::build_R() {
        double tol = 1.0e-20;
        size_t npts = get_npts();
        arma::vec temp(npts);
        arma::vec pi(get_pi());
        arma::vec rho(get_rho());
        #pragma parallel for schedule(static) \
	                     private(p) \
                             shared(npts, pi, rho, temp)
           for (int p = 0; p < npts; p++) {
               temp(p) = 4.0 * pi(p) / ( rho(p) * rho(p) );
           }
        set_R(temp);
   }

   void MCPDFT::translate() {
     translate_density();
     if ( is_gga_ ) {
        translate_density_gradients();
     }
   }

   void MCPDFT::translate_density() {
      double tol = 1.0e-20;
      size_t npts = get_npts();
      arma::vec rho_vec(get_rho());
      arma::vec pi_vec(get_pi());
      arma::vec R_vec(get_R());
      arma::vec W(get_w());
      arma::vec tr_rhoa(npts);
      arma::vec tr_rhob(npts);
      double dum_a = 0.0;
      double dum_b = 0.0;
      double dum_tot = 0.0;
      double rho = 0.0;
      double pi = 0.0;
      double zeta = 0.0;
      double R = 0.0;
      #pragma parallel for schedule(static) \
                           default(shared) \
                           private(p, zeta, pi, rho) \
                           reduction(+:dum_a, dum_b, dum_tot)
         for (int p = 0; p < npts; p++) {
             rho = rho_vec(p);
             pi = pi_vec(p);
             zeta = 0.0;
             R = 0.0;
             if ( !(rho < tol) && !(pi < 0.0) ) {
                R = R_vec(p);
                if ( (1.0 - R) > tol ) {
                   zeta = sqrt(1.0 - R);
                }else{
                     zeta = 0.0;
                }
                tr_rhoa(p) = (1.0 + zeta) * (rho/2.0);
                tr_rhob(p) = (1.0 - zeta) * (rho/2.0);
             }else {
                    tr_rhoa(p) = 0.0;
                    tr_rhob(p) = 0.0;
             }
             dum_a += tr_rhoa(p) * W(p);
             dum_b += tr_rhob(p) * W(p);
             dum_tot += ( tr_rhoa(p) + tr_rhob(p) ) * W(p) ;
         }
      set_tr_rhoa(tr_rhoa);
      set_tr_rhob(tr_rhob);

      printf("\n");
      printf("  Integrated translated total density = %20.12lf\n",dum_tot);
      printf("  Integrated translated alpha density = %20.12lf\n",dum_a);
      printf("  Integrated translated beta density  = %20.12lf\n",dum_b);
      printf("\n");
   }
 
   void MCPDFT::translate_density_gradients() {
     double tol = 1.0e-20;
     size_t npts = get_npts();
     arma::vec rho_vec(get_rho());
     arma::vec pi_vec(get_pi());
     arma::vec R_vec(get_R());
     double rho = 0.0;
     double pi = 0.0;
     double zeta = 0.0;
     double R = 0.0;
     double rho_x = 0.0;
     double rho_y = 0.0;
     double rho_z = 0.0;
     arma::vec rho_a_x(get_rhoa_x());
     arma::vec rho_a_y(get_rhoa_y());
     arma::vec rho_a_z(get_rhoa_z());
     arma::vec rho_b_x(get_rhob_x());
     arma::vec rho_b_y(get_rhob_y());
     arma::vec rho_b_z(get_rhob_z());
     arma::vec tr_rho_a_x(npts,  arma::fill::zeros);
     arma::vec tr_rho_b_x(npts,  arma::fill::zeros);
     arma::vec tr_rho_a_y(npts,  arma::fill::zeros);
     arma::vec tr_rho_b_y(npts,  arma::fill::zeros);
     arma::vec tr_rho_a_z(npts,  arma::fill::zeros);
     arma::vec tr_rho_b_z(npts,  arma::fill::zeros);
     arma::vec tr_sigma_aa(npts, arma::fill::zeros);
     arma::vec tr_sigma_ab(npts, arma::fill::zeros);
     arma::vec tr_sigma_bb(npts, arma::fill::zeros);
     for (int p = 0; p < npts; p++) {
         rho = rho_vec(p);
         pi = pi_vec(p);
         rho_x = rho_a_x(p) + rho_b_x(p);
         rho_y = rho_a_y(p) + rho_b_y(p);
         rho_z = rho_a_z(p) + rho_b_z(p);
         zeta = 0.0;
         R = 0.0;
         if ( !(rho < tol) && !(pi < 0.0) ) {
            R = R_vec(p);
            if ( (1.0 - R) > tol )  {
               zeta = sqrt(1.0 - R);
            }else{
                 zeta = 0.0;
            }
            tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0);
            tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0);

            tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0);
            tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0);

            tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0);
            tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0);
         }else {
               tr_rho_a_x(p) = 0.0;
               tr_rho_b_x(p) = 0.0;
               tr_rho_a_y(p) = 0.0;
               tr_rho_b_y(p) = 0.0;
               tr_rho_a_z(p) = 0.0;
               tr_rho_b_z(p) = 0.0;
         }
         tr_sigma_aa(p) = (tr_rho_a_x(p) * tr_rho_a_x(p)) + (tr_rho_a_y(p) * tr_rho_a_y(p)) + (tr_rho_a_z(p) * tr_rho_a_z(p));
         tr_sigma_ab(p) = (tr_rho_a_x(p) * tr_rho_b_x(p)) + (tr_rho_a_y(p) * tr_rho_b_y(p)) + (tr_rho_a_z(p) * tr_rho_b_z(p));
         tr_sigma_bb(p) = (tr_rho_b_x(p) * tr_rho_b_x(p)) + (tr_rho_b_y(p) * tr_rho_b_y(p)) + (tr_rho_b_z(p) * tr_rho_b_z(p));
     }
     set_tr_sigma_aa(tr_sigma_aa);
     set_tr_sigma_ab(tr_sigma_ab);
     set_tr_sigma_bb(tr_sigma_bb);
   }

   void MCPDFT::fully_translate(){
        double tol = 1.0e-20;
        size_t npts = get_npts();

        arma::vec rho_vec(get_rho());
        arma::vec pi_vec(get_pi());
        arma::vec R_vec(get_R());
        arma::vec W(get_w());
        arma::vec tr_rhoa(npts);
        arma::vec tr_rhob(npts);
        arma::vec tr_rho(npts);

        double const R0 = 0.9;
        double const R1 = 1.15;
        double const A = -475.60656009;
        double const B = -379.47331922;
        double const C = -85.38149682;

        double temp_tot = 0.0;
        double temp_a = 0.0;
        double temp_b = 0.0;
        for (int p = 0; p < npts; p++) {
            double zeta = 0.0;
            double R = 0.0;
            double rho = rho_vec(p);
            double pi = pi_vec(p);
            double DelR = R_vec(p) - R1;
            if ( !(rho < tol) && !(pi < 0.0) ) {
               R = R_vec(p);
               if ( ((1.0 - R) > tol) && ( R < R0 ) ) {
                  zeta = sqrt(1.0 - R);
               }else if( !(R < R0) && !(R > R1) ) {
                       zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);
               }else if( R > R1 ) {
                       zeta = 0.0;
               }
               tr_rhoa(p) = (1.0 + zeta) * (rho/2.0);
               tr_rhob(p) = (1.0 - zeta) * (rho/2.0);
            }else{
                 tr_rhoa(p) = 0.0;
                 tr_rhob(p) = 0.0;
            }
            temp_a += tr_rhoa(p) * W(p);
            temp_b += tr_rhob(p) * W(p);
            temp_tot += ( tr_rhob(p) + tr_rhoa(p) ) * W(p);
        }
        set_tr_rhoa(tr_rhoa);
        set_tr_rhob(tr_rhob);

        printf("\n");
        printf("      Integrated fully translated total density = %20.12lf\n",temp_tot);
        printf("      Integrated fully translated alpha density = %20.12lf\n",temp_a);
        printf("      Integrated fully translated beta density  = %20.12lf\n",temp_b);
        printf("\n");

        //if ( is_gga_ || is_meta_ ) {
        // arma::vec tr_rho_a_x(npts);
        // arma::vec tr_rho_b_x(npts);

        // arma::vec tr_rho_a_y(npts);
        // arma::vec tr_rho_b_y(npts);

        // arma::vec tr_rho_a_z(npts);
        // arma::vec tr_rho_b_z(npts);

        // arma::vec tr_sigma_aa(npts);
        // arma::vec tr_sigma_bb(npts);
        // arma::vec tr_sigma_ab(npts);

        // arma::mat temp_tr_rho1_a(npts,3);
        // arma::mat temp_tr_rho1_b(npts,3);
        // for (int p = 0; p < npts; p++) {

        //     double rho_x = rho_a_x_(p) + rho_b_x_(p);
        //     double rho_y = rho_a_y_(p) + rho_b_y_(p);
        //     double rho_z = rho_a_z_(p) + rho_b_z_(p);

        //     double rho = rho_p(p);
        //     double pi = pi_p(p);
        //     double DelR = R_p(p) - R1;

        //     double zeta = 0.0;
        //     double R = 0.0;

        //     // if ( !(rho < tol) && !(pi < tol) ) {
        //     if ( !(rho < tol) && !(pi < 0.0) ) {

        //         R = R_p(p);
        //         // R = tanh(R);

        //         if ( ((1.0 - R) > tol) && ( R < R0 ) ) {

        //            zeta = sqrt(1.0 - R);

        //            tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0) + (R * rho_x) / (2.0*zeta) - pi_x_(p) / (rho*zeta);
        //            tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0) - (R * rho_x) / (2.0*zeta) + pi_x_(p) / (rho*zeta);

        //            tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0) + (R * rho_y) / (2.0*zeta) - pi_y_(p) / (rho*zeta);
        //            tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0) - (R * rho_y) / (2.0*zeta) + pi_y_(p) / (rho*zeta);

        //            tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0) + (R * rho_z) / (2.0*zeta) - pi_z_(p) / (rho*zeta);
        //            tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0) - (R * rho_z) / (2.0*zeta) + pi_z_(p) / (rho*zeta);

        //         }else if( !(R < R0) && !(R > R1) ) {

        //                 zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

        //                 tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0)
        //                                + (A * pow(DelR, 4.0)) * ( (10.0 * pi_x_(p) / rho) - (5.0 * R * rho_x) )
        //                                + (B * pow(DelR, 3.0)) * ( (8.0  * pi_x_(p) / rho) - (4.0 * R * rho_x) )
        //                                + (C * pow(DelR, 2.0)) * ( (6.0  * pi_x_(p) / rho) - (3.0 * R * rho_x) );

        //                 tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0)
        //                                + (A * pow(DelR, 4.0)) * (-(10.0 * pi_x_(p) / rho) + (5.0 * R * rho_x) )
        //                                + (B * pow(DelR, 3.0)) * (-(8.0  * pi_x_(p) / rho) + (4.0 * R * rho_x) )
        //                                + (C * pow(DelR, 2.0)) * (-(6.0  * pi_x_(p) / rho) + (3.0 * R * rho_x) );

        //                 tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0)
        //                                + (A * pow(DelR, 4.0)) * ( (10.0 * pi_y_(p) / rho) - (5.0 * R * rho_y) )
        //                                + (B * pow(DelR, 3.0)) * ( (8.0  * pi_y_(p) / rho) - (4.0 * R * rho_y) )
        //                                + (C * pow(DelR, 2.0)) * ( (6.0  * pi_y_(p) / rho) - (3.0 * R * rho_y) );

        //                 tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0)
        //                                + (A * pow(DelR, 4.0)) * (-(10.0 * pi_y_(p) / rho) + (5.0 * R * rho_y) )
        //                                + (B * pow(DelR, 3.0)) * (-(8.0  * pi_y_(p) / rho) + (4.0 * R * rho_y) )
        //                                + (C * pow(DelR, 2.0)) * (-(6.0  * pi_y_(p) / rho) + (3.0 * R * rho_y) );

        //                 tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0)
        //                                + (A * pow(DelR, 4.0)) * ( (10.0 * pi_z_(p) / rho) - (5.0 * R * rho_z) )
        //                                + (B * pow(DelR, 3.0)) * ( (8.0  * pi_z_(p) / rho) - (4.0 * R * rho_z) )
        //                                + (C * pow(DelR, 2.0)) * ( (6.0  * pi_z_(p) / rho) - (3.0 * R * rho_z) );

        //                 tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0)
        //                                 + (A * pow(DelR, 4.0)) * (-(10.0 * pi_z_(p) / rho) + (5.0 * R * rho_z) )
        //                                 + (B * pow(DelR, 3.0)) * (-(8.0  * pi_z_(p) / rho) + (4.0 * R * rho_z) )
        //                                 + (C * pow(DelR, 2.0)) * (-(6.0  * pi_z_(p) / rho) + (3.0 * R * rho_z) );
        //         }else if( R > R1 ) {
        //                 zeta = 0.0;

        //                 tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0);
        //                 tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0);

        //                 tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0);
        //                 tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0);

        //                 tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0);
        //                 tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0);
        //         }

        //     }else{

        //         tr_rho_a_x(p) = 0.0;
        //         tr_rho_b_x(p) = 0.0;

        //         tr_rho_a_y(p) = 0.0;
        //         tr_rho_b_y(p) = 0.0;

        //         tr_rho_a_z(p) = 0.0;
        //         tr_rho_b_z(p) = 0.0;
        //     }

        //     temp_tr_rho1_a(p,0) = tr_rho_a_x(p); 
        //     temp_tr_rho1_a(p,1) = tr_rho_a_y(p); 
        //     temp_tr_rho1_a(p,2) = tr_rho_a_z(p); 
        //     temp_tr_rho1_b(p,0) = tr_rho_b_x(p); 
        //     temp_tr_rho1_b(p,1) = tr_rho_b_y(p); 
        //     temp_tr_rho1_b(p,2) = tr_rho_b_z(p); 

        //     tr_sigma_aa(p) = (tr_rho_a_x(p) * tr_rho_a_x(p)) + (tr_rho_a_y(p) * tr_rho_a_y(p)) + (tr_rho_a_z(p) * tr_rho_a_z(p));
        //     tr_sigma_ab(p) = (tr_rho_a_x(p) * tr_rho_b_x(p)) + (tr_rho_a_y(p) * tr_rho_b_y(p)) + (tr_rho_a_z(p) * tr_rho_b_z(p));
        //     tr_sigma_bb(p) = (tr_rho_b_x(p) * tr_rho_b_x(p)) + (tr_rho_b_y(p) * tr_rho_b_y(p)) + (tr_rho_b_z(p) * tr_rho_b_z(p));

        //  //}

        // }
        // set_tr_rho1_a(temp_tr_rho1_a);
        // set_tr_rho1_b(temp_tr_rho1_b);
   }
   void MCPDFT::build_opdm() {
      // fetching the number of basis functions
      int nbfs;
      nbfs = get_nbfs();
 
      // getting the AO->MO transformation matrix C
      arma::mat ca(get_cmat());
      arma::mat cb(get_cmat());
 
      // building the 1-electron reduced density matrices (1RDMs)
      arma::mat D1a(nbfs, nbfs, arma::fill::zeros);
      arma::mat D1b(nbfs, nbfs, arma::fill::zeros);
      for (int mu = 0; mu < nbfs; mu++) { 
          for (int nu = 0; nu < nbfs; nu++) { 
              double duma = 0.0;
              double dumb = 0.0;
              for (int i = 0; i < nbfs/2; i++) { 
                   duma += ca(mu, i) * ca(nu, i);
                   dumb += cb(mu, i) * cb(nu, i);
              }
              D1a(mu, nu) = duma;
              D1b(mu, nu) = dumb;
         }
      }
      // D1a(0,0) = 1.0;
      // D1b(0,0) = 1.0;
      D1a.print("D1a = ");
      D1b.print("D1b = ");
      set_D1a(D1a);
      set_D1b(D1b);
   }
 
   void MCPDFT::build_tpdm() {
      // fetching the number of basis functions
      int nbfs  = get_nbfs();
      int nbfs2 = nbfs * nbfs;
 
      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());
 
      arma::mat D2ab(nbfs2, nbfs2, arma::fill::zeros);
      D2ab = arma::kron(D1a,D1b);
      // D2ab.print("D2ab = ");
      set_D2ab(D2ab);
   }
 
   bool MCPDFT::is_gga() const { return is_gga_; }
   size_t MCPDFT::get_npts() const { return npts_; }
   int    MCPDFT::get_nbfs() const { return nbfs_; }
   arma::vec MCPDFT::get_w() const { return w_; }
   arma::vec MCPDFT::get_x() const { return x_; }
   arma::vec MCPDFT::get_y() const { return y_; }
   arma::vec MCPDFT::get_z() const { return z_; }
   arma::mat MCPDFT::get_phi() const { return phi_; }
   arma::mat MCPDFT::get_phi_x() const { return phi_x_; }
   arma::mat MCPDFT::get_phi_y() const { return phi_y_; }
   arma::mat MCPDFT::get_phi_z() const { return phi_z_; }
   double MCPDFT::get_eref() const { return eref_; }
   double MCPDFT::get_eclass()  const { return eclass_; }
   arma::mat MCPDFT::get_cmat() const { return cmat_; }
   arma::mat MCPDFT::get_D1a()  const { return D1a_ ; }
   arma::mat MCPDFT::get_D1b()  const { return D1b_ ; }
   arma::mat MCPDFT::get_D2ab()  const { return D2ab_ ; }
   arma::vec MCPDFT::get_rhoa() const { return rho_a_; }
   arma::vec MCPDFT::get_rhoa_x() const { return rho_a_x_; }
   arma::vec MCPDFT::get_rhoa_y() const { return rho_a_y_; }
   arma::vec MCPDFT::get_rhoa_z() const { return rho_a_z_; }
   arma::vec MCPDFT::get_rhob() const { return rho_b_; }
   arma::vec MCPDFT::get_rhob_x() const { return rho_b_x_; }
   arma::vec MCPDFT::get_rhob_y() const { return rho_b_y_; }
   arma::vec MCPDFT::get_rhob_z() const { return rho_b_z_; }
   arma::vec MCPDFT::get_rho() const { return rho_; }
   arma::vec MCPDFT::get_tr_rhoa() const { return tr_rho_a_; }
   arma::vec MCPDFT::get_tr_rhob() const { return tr_rho_b_; }
   arma::vec MCPDFT::get_tr_rho() const { return tr_rho_; }
   arma::vec MCPDFT::get_pi() const { return pi_; }
   arma::vec MCPDFT::get_R() const { return R_; }
   arma::vec MCPDFT::get_sigma_aa() const { return sigma_aa_; }
   arma::vec MCPDFT::get_sigma_ab() const { return sigma_ab_; }
   arma::vec MCPDFT::get_sigma_bb() const { return sigma_bb_; }
   arma::vec MCPDFT::get_tr_sigma_aa() const { return tr_sigma_aa_; }
   arma::vec MCPDFT::get_tr_sigma_ab() const { return tr_sigma_ab_; }
   arma::vec MCPDFT::get_tr_sigma_bb() const { return tr_sigma_bb_; }

   void MCPDFT::set_npts(const size_t npts) { npts_ = npts; }
   void MCPDFT::set_nbfs(const int nbfs)    { nbfs_ = nbfs; }
   void MCPDFT::set_w(const arma::vec &w) { w_ = w; }
   void MCPDFT::set_x(const arma::vec &x) { x_ = x; }
   void MCPDFT::set_y(const arma::vec &y) { y_ = y; }
   void MCPDFT::set_z(const arma::vec &z) { z_ = z; }
   void MCPDFT::set_phi(const arma::mat &phi) { phi_ = phi; }
   void MCPDFT::set_phi_x(const arma::mat &phi_x) { phi_x_ = phi_x; }
   void MCPDFT::set_phi_y(const arma::mat &phi_y) { phi_y_ = phi_y; }
   void MCPDFT::set_phi_z(const arma::mat &phi_z) { phi_z_ = phi_z; }
   void MCPDFT::set_eref(const double eref) { eref_ = eref; }
   void MCPDFT::set_eclass(const double eclass) { eclass_ = eclass; }
   void MCPDFT::set_cmat(const arma::mat &cmat) { cmat_ = cmat; }
   void MCPDFT::set_D1a(const arma::mat &D1a) { D1a_ = D1a; }
   void MCPDFT::set_D1b(const arma::mat &D1b) { D1b_ = D1b; }
   void MCPDFT::set_D2ab(const arma::mat &D2ab) { D2ab_ = D2ab; }
   void MCPDFT::set_rhoa(const arma::vec &rhoa) { rho_a_ = rhoa; }
   void MCPDFT::set_rhoa_x(const arma::vec &rhoa_x) { rho_a_x_ = rhoa_x; }
   void MCPDFT::set_rhoa_y(const arma::vec &rhoa_y) { rho_a_y_ = rhoa_y; }
   void MCPDFT::set_rhoa_z(const arma::vec &rhoa_z) { rho_a_z_ = rhoa_z; }
   void MCPDFT::set_rhob(const arma::vec &rhob) { rho_b_ = rhob; }
   void MCPDFT::set_rhob_x(const arma::vec &rhob_x) { rho_b_x_ = rhob_x; }
   void MCPDFT::set_rhob_y(const arma::vec &rhob_y) { rho_b_y_ = rhob_y; }
   void MCPDFT::set_rhob_z(const arma::vec &rhob_z) { rho_b_z_ = rhob_z; }
   void MCPDFT::set_rho(const arma::vec &rho) { rho_ = rho; }
   void MCPDFT::set_tr_rhoa(const arma::vec &tr_rhoa) { tr_rho_a_ = tr_rhoa; }
   void MCPDFT::set_tr_rhoa_x(const arma::vec &tr_rhoa_x) { tr_rho_a_x_ = tr_rhoa_x; }
   void MCPDFT::set_tr_rhoa_y(const arma::vec &tr_rhoa_y) { tr_rho_a_y_ = tr_rhoa_y; }
   void MCPDFT::set_tr_rhoa_z(const arma::vec &tr_rhoa_z) { tr_rho_a_z_ = tr_rhoa_z; }
   void MCPDFT::set_tr_rhob(const arma::vec &tr_rhob) { tr_rho_b_ = tr_rhob; }
   void MCPDFT::set_tr_rhob_x(const arma::vec &tr_rhob_x) { tr_rho_b_x_ = tr_rhob_x; }
   void MCPDFT::set_tr_rhob_y(const arma::vec &tr_rhob_y) { tr_rho_b_y_ = tr_rhob_y; }
   void MCPDFT::set_tr_rhob_z(const arma::vec &tr_rhob_z) { tr_rho_b_z_ = tr_rhob_z; }
   void MCPDFT::set_tr_rho(const arma::vec &tr_rho) { tr_rho_ = tr_rho; }
   void MCPDFT::set_pi(const arma::vec &pi) { pi_ = pi; }
   void MCPDFT::set_pi_x(const arma::vec &pi_x) { pi_x_ = pi_x; }
   void MCPDFT::set_pi_y(const arma::vec &pi_y) { pi_y_ = pi_y; }
   void MCPDFT::set_pi_z(const arma::vec &pi_z) { pi_z_ = pi_z; }
   void MCPDFT::set_R(const arma::vec &R) { R_ = R; }
   void MCPDFT::set_sigma_aa(const arma::vec &sigma_aa) { sigma_aa_ = sigma_aa; }
   void MCPDFT::set_sigma_ab(const arma::vec &sigma_ab) { sigma_ab_ = sigma_ab; }
   void MCPDFT::set_sigma_bb(const arma::vec &sigma_bb) { sigma_bb_ = sigma_bb; }
   void MCPDFT::set_tr_sigma_aa(const arma::vec &tr_sigma_aa) { tr_sigma_aa_ = tr_sigma_aa; }
   void MCPDFT::set_tr_sigma_ab(const arma::vec &tr_sigma_ab) { tr_sigma_ab_ = tr_sigma_ab; }
   void MCPDFT::set_tr_sigma_bb(const arma::vec &tr_sigma_bb) { tr_sigma_bb_ = tr_sigma_bb; }

   void MCPDFT::print_banner() const {
      printf("\n******************************************************************\n");
      printf("*                                                                *\n");
      printf("*                           OpenRDM:                             *\n");
      printf("*                                                                *\n");
      printf("*                 An open-source library for                     *\n");
      printf("*    reduced-density matrix-based analysis and computation       *\n");
      printf("*                                                                *\n");
      printf("*                     Mohammad Mostafanejad                      *\n");
      printf("*                   Florida State University                     *\n");
      printf("*                                                                *\n");
      printf("******************************************************************\n");


      printf("\n           Please cite the following article(s):\n\n");

      printf("    # M. Mostafanejad and A. E. DePrince III\n");
      printf("      J. Chem. Theory Comput. 15, 290-302 (2019).\n");
   }
}
