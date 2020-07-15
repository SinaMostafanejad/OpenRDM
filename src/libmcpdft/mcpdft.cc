#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string>
#include "mcpdft.h"
#include "openrdmConfig.h"

#include "HDF5Client.h"
#include "HDF5Utility.h"
#include "TOC.h"

#ifdef WITH_OPENMP // _OPENMP
   #include <omp.h>
#endif

namespace mcpdft {

   MCPDFT::MCPDFT(std::string test_case)  { common_init(test_case); }
   MCPDFT::MCPDFT() { common_init(); }
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

   void MCPDFT::common_init() {
      print_banner();
      /* TODO: make some setter and getter for this to work at the moment
         until the parser gets implemented. */
      is_ao_ = false;
      is_gga_ = true;
      is_sparse_ = true;
      is_active_ = false;

      HDF5Utility* h5utl = new HDF5Utility();

      /* reading naos, nmos and npts from data.h5 HDF5 file */	   
      size_t nao{0}, nmo{0}, npts{0}, nnz{0};
      h5utl->read_basics(nao, nmo, npts);
//      std::printf("NAO, NMO, NPTS = %ld, %ld, %ld\n",nao,nmo,npts);
      set_nao(nao);
      set_nmo(nmo);
      set_npts(npts);

      /* reading active space details from data.h5 HDF5 file */
      size_t nactele{0}, nactorb{0};
      size_t ncore{0}, nfrz{0}, nocc{0}, nvir{0};
      h5utl->read_active_space_vars(nactele,
		                    nactorb,
				    ncore,
				    nfrz,
				    nocc,
				    nvir);

//      std::printf("NCASELE, NCASOBR, NCORE, NFRZ, NOCC, NVIR = %d , %d, %d, %d, %d, %d\n",
//                  nactele,nactorb,ncore,nfrz,nocc,nvir);


      /* reading AO to MO transformation matrix from data.h5 HDF5 file */
      arma::mat Cmat(nmo, nao, arma::fill::zeros);
      h5utl->read_AO2MO_Cmat(Cmat);
      Cmat = Cmat.t();
//      Cmat.print();
      set_cmat(Cmat);

      /* reading classical nuclear repulsion energy from data.h5 HDF5 file */
      double tmp_enuc{0.0};
      h5utl->read_enuc(tmp_enuc);
      set_enuc(tmp_enuc);

      /* reading core Hamiltonian matrix in AO/MO basis from data.h5 HDF5 file */
      arma::mat Hcore(nmo, nmo, arma::fill::zeros);
      h5utl->read_hcore(Hcore, is_ao_);
      Hcore = Hcore.t();
//      Hcore.print();
      set_hcore(Hcore);

      /* reading (Coulomb) Hartree interaction matrices in AO/MO basis from data.h5 HDF5 file */
      arma::mat Ja(nmo, nmo, arma::fill::zeros);
      arma::mat Jb(nmo, nmo, arma::fill::zeros);
      h5utl->read_Hartree_Jmats(Ja, Jb, is_ao_);
      Ja = Ja.t();
      Jb = Jb.t();
//      Ja.print();
//      Jb.print();
      set_ja(Ja);
      set_jb(Jb);

      /* reading Cartesian grids and weights from data.h5 HDF5 file */
      arma::vec W(npts, arma::fill::zeros);
      arma::vec X(npts, arma::fill::zeros);
      arma::vec Y(npts, arma::fill::zeros);
      arma::vec Z(npts, arma::fill::zeros);

      h5utl->read_grids(W, X, Y, Z);
//      W.print("W = ");
//      X.print("X = ");
//      Y.print("Y = ");
//      Z.print("Z = ");
      set_w(W);
      set_x(X);
      set_y(Y);
      set_z(Z);

      /* reading AOs/MOs matrices calculated on grid points from data.h5 HDF5 file */
      // fixing the is_ao to false for now.
      arma::mat phi(nmo, npts, arma::fill::zeros);
      arma::mat phi_x(nmo, npts,  arma::fill::zeros);
      arma::mat phi_y(nmo, npts,  arma::fill::zeros);
      arma::mat phi_z(nmo, npts,  arma::fill::zeros);
      h5utl->read_superphi(phi,
		           phi_x,
			   phi_y,
			   phi_z,
			   is_ao_);
//      phi.t().print("PHI = ");
//      phi_x.t().print("PHI_X = ");
//      phi_y.t().print("PHI_Y = ");
//      phi_z.t().print("PHI_Z = ");
      set_phi(phi);
      set_phi_x(phi_x);
      set_phi_y(phi_y);
      set_phi_z(phi_z);
//      set_phi(phi.t());
//      set_phi_x(phi_x.t());
//      set_phi_y(phi_y.t());
//      set_phi_z(phi_z.t());

      arma::mat d1a(nmo, nmo, arma::fill::zeros);
      arma::mat d1b(nmo, nmo, arma::fill::zeros);
      arma::mat d2ab(nmo*nmo, nmo*nmo, arma::fill::zeros);

      HDF5Client* h5client = new HDF5Client();
      HDF5Client::factory_mode mode(HDF5Client::factory_mode::READ);
      h5client->factory_client(H5D_CONTIGUOUS,mode,d1a,d1b,d2ab);
      delete h5client;

//      d1a.print("D1a = ");
//      d1b.print("D1b = ");
//      d2ab.t().print("D2ab = ");

      set_D1a(d1a);
      set_D1b(d1b);
      set_D2ab(d2ab);

      delete h5utl;
   }

   void MCPDFT::build_rho() {
      const bool is_sparse = MCPDFT::is_sparse();

      if(!is_sparse) { /* using dens RDMs*/
         build_density_functions();
         if ( is_gga_ ) {
            build_density_gradients();
         }
      }else{ /* using sparse RDMs*/
         build_sparse_density_functions();
         if ( is_gga_ ) {
            build_sparse_density_gradients();
         }
      }
   }

   void MCPDFT::build_density_functions() {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
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
      double tmp_d1a{0.0}, tmp_d1b{0.0};
      #ifdef WITH_OPENMP
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
         printf("                   *** Warning ***\n");
         printf("   Calculating the density (gradients) using OpenMP");
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      #endif

      #pragma omp parallel default(shared) \
                           private(p, mu, nu,\
       			    tmp_d1a, tmp_d1b)
      {  
         for(int mu = 0; mu < nbfs; mu++) {
            for(int nu = 0; nu < nbfs; nu++) {
               tmp_d1a = D1a(mu, nu);
               tmp_d1b = D1b(mu, nu);
               if (fabs(tmp_d1a) < tol && fabs(tmp_d1b) < tol )
                  continue;
               #pragma omp for schedule(static)
                  for(int p = 0; p < npts; p++) {
                     rhoa(p) += tmp_d1a * phi(mu, p) * phi(nu, p);
                     rhob(p) += tmp_d1b * phi(mu, p) * phi(nu, p);
                  }
            }
         }
      } /* end of omp parallel region */
      set_rhoa(rhoa);
      set_rhob(rhob);

      for(p = 0; p < npts; p++) {
         rho(p) = rhoa(p) + rhob(p);
      } 
      set_rho(rho);

      dum_a = arma::dot(rhoa, W);
      dum_b = arma::dot(rhob, W);
      dum_tot = arma::dot(rho, W);
      printf("\n");
      printf("  Integrated alpha density = %20.12lf\n",dum_a);
      printf("  Integrated beta density  = %20.12lf\n",dum_b);
      printf("  Integrated total density = %20.12lf\n",dum_tot);
      printf("\n");
   }

   void MCPDFT::build_sparse_density_functions() {
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::vec W(get_w());
      arma::vec rhoa(npts, arma::fill::zeros);
      arma::vec rhob(npts, arma::fill::zeros);
      arma::vec rho(npts, arma::fill::zeros);
      bool is_active = MCPDFT::is_active();

      size_t nnz{0};
      HDF5Utility* h5utl = new HDF5Utility();

      h5utl->read_nnz(nnz, is_active, "D1A");
      // std::printf("nnz = %6.lu",nnz);

      arma::vec val_a(nnz, arma::fill::zeros);
      arma::vec val_b(nnz, arma::fill::zeros);
      arma::Col<int> row_idx_a(nnz, arma::fill::zeros);
      arma::Col<int> row_idx_b(nnz, arma::fill::zeros);
      arma::Col<int> col_idx_a(nnz, arma::fill::zeros);
      arma::Col<int> col_idx_b(nnz, arma::fill::zeros);

      h5utl->read_sparse_coo_opdm(val_a, row_idx_a, col_idx_a, is_active, "D1A");
      h5utl->read_sparse_coo_opdm(val_b, row_idx_b, col_idx_b, is_active, "D1B");

      delete h5utl;

      for (size_t n = 0; n < nnz; n++) {
	 size_t i_a = row_idx_a(n);
	 size_t i_b = row_idx_b(n);
	 size_t j_a = col_idx_a(n);
	 size_t j_b = col_idx_b(n);
	 /* since we stored upper triangular part of the symmetric 1-RDM,
	    we need a factor of 2.0 for off-diagonal elements here */
	 double s_a = (i_a != j_a) ? 2.0 : 1.0;
	 double s_b = (i_b != j_b) ? 2.0 : 1.0;
	 double tmp_val_a = val_a(n);
	 double tmp_val_b = val_b(n);
         for(size_t p = 0; p < npts; p++) {
	    rhoa(p) += s_a * tmp_val_a * phi(i_a, p) * phi(j_a, p);
	    rhob(p) += s_b * tmp_val_b * phi(i_b, p) * phi(j_b, p);
	 }
      }
      set_rhoa(rhoa);
      set_rhob(rhob);

      for(size_t p = 0; p < npts; p++) {
         rho(p) = rhoa(p) + rhob(p);
      } 
      set_rho(rho);

      double dum_a = arma::dot(rhoa, W);
      double dum_b = arma::dot(rhob, W);
      double dum_tot = arma::dot(rho, W);
      printf("\n");
      printf("  Integrated alpha density (sparse) = %20.12lf\n",dum_a);
      printf("  Integrated beta density  (sparse) = %20.12lf\n",dum_b);
      printf("  Integrated total density (sparse) = %20.12lf\n",dum_tot);
      printf("\n");
   } 
#if 0
   void MCPDFT::build_density_functions() {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
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
      double tmp_d1a{0.0}, tmp_d1b{0.0};

      #ifdef WITH_OPENMP
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
         printf("                   *** Warning ***\n");
         printf("   Calculating the density (gradients) using OpenMP");
         printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      #endif

      #pragma omp parallel default(shared) \
                           private(p, mu, nu, \
			           tmp_d1a, tmp_d1b)
      {  
         for(int mu = 0; mu < nbfs; mu++) {
            for(int nu = 0; nu < nbfs; nu++) {
	       tmp_d1a = D1a(mu, nu);
	       tmp_d1b = D1b(mu, nu);
	       if (fabs(tmp_d1a) < tol && fabs(tmp_d1b) < tol )
	          continue;
               #pragma omp for schedule(static)
                  for(int p = 0; p < npts; p++) {
                     rhoa(p) += tmp_d1a * phi(mu, p) * phi(nu, p);
                     rhob(p) += tmp_d1b * phi(mu, p) * phi(nu, p);
                     rho(p) += rhoa(p) + rhob(p);
	          }
            }
         }

      } /* end of omp parallel region */
      set_rhoa(rhoa);
      set_rhob(rhob);
      set_rho(rho);
         //#pragma omp for schedule(static) \
         //                reduction(+:dum_a, dum_b, dum_tot) \
	 //                nowait
         //   for(int p = 0; p < npts; p++) {
         //         dum_a += rhoa(p) * W(p);
         //         dum_b += rhob(p) * W(p);
         //         dum_tot += ( rhoa(p) + rhob(p) ) * W(p) ;
         //   } /* end of omp parallel for loop */
      dum_a = arma::dot(rhoa, W);
      dum_b = arma::dot(rhob, W);
      dum_tot = dum_a + dum_b;
      printf("\n");
      printf("  Integrated total density = %20.12lf\n",dum_tot);
      printf("  Integrated alpha density = %20.12lf\n",dum_a);
      printf("  Integrated beta density  = %20.12lf\n",dum_b);
      printf("\n");
   }
#endif

   void MCPDFT::build_density_gradients() {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());
      int p{0}, nu{0}, sigma{0};
      double tmp_d1a{0.0}, tmp_d1b{0.0};
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
                           private(p, nu, sigma,\
		                   tmp_d1a, tmp_d1b)
      {
                    for (int sigma = 0; sigma < nbfs; sigma++) {
                        for (int nu = 0; nu < nbfs; nu++) {
	                   tmp_d1a = D1a(sigma, nu);
	                   tmp_d1b = D1b(sigma, nu);
	                   if (fabs(tmp_d1a) < tol && fabs(tmp_d1b) < tol )
	                      continue;
                           #pragma omp for schedule(static) 
                              for (int p = 0; p < npts; p++) {
                                 rho_a_x(p) += ( phi_x(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_x(nu, p) ) * tmp_d1a;
                                 rho_b_x(p) += ( phi_x(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_x(nu, p) ) * tmp_d1b;
                                 rho_a_y(p) += ( phi_y(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_y(nu, p) ) * tmp_d1a;
                                 rho_b_y(p) += ( phi_y(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_y(nu, p) ) * tmp_d1b;
                                 rho_a_z(p) += ( phi_z(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_z(nu, p) ) * tmp_d1a;
                                 rho_b_z(p) += ( phi_z(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_z(nu, p) ) * tmp_d1b;
			      }
                        }
                    }
      }
      set_rhoa_x(rho_a_x);
      set_rhob_x(rho_b_x);
      set_rhoa_y(rho_a_y);
      set_rhob_y(rho_b_y);
      set_rhoa_z(rho_a_z);
      set_rhob_z(rho_b_z);
      // rho_a_x.print();

      for (int p = 0; p < npts; p++) {
         sigma_aa(p) = ( rho_a_x(p) * rho_a_x(p) ) +  ( rho_a_y(p) * rho_a_y(p) ) + ( rho_a_z(p) * rho_a_z(p) );
         sigma_bb(p) = ( rho_b_x(p) * rho_b_x(p) ) +  ( rho_b_y(p) * rho_b_y(p) ) + ( rho_b_z(p) * rho_b_z(p) );
         sigma_ab(p) = ( rho_a_x(p) * rho_b_x(p) ) +  ( rho_a_y(p) * rho_b_y(p) ) + ( rho_a_z(p) * rho_b_z(p) );
      }
      set_sigma_aa(sigma_aa);
      set_sigma_ab(sigma_ab);
      set_sigma_bb(sigma_bb);
   }

   void MCPDFT::build_sparse_density_gradients() {
      size_t npts = get_npts();
      arma::mat phi(get_phi());
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
      bool is_active = MCPDFT::is_active();

      size_t nnz{0};
      HDF5Utility* h5utl = new HDF5Utility();

      h5utl->read_nnz(nnz, is_active, "D1A");
      // std::printf("nnz = %6.lu",nnz);

      arma::vec val_a(nnz, arma::fill::zeros);
      arma::vec val_b(nnz, arma::fill::zeros);
      arma::Col<int> row_idx_a(nnz, arma::fill::zeros);
      arma::Col<int> row_idx_b(nnz, arma::fill::zeros);
      arma::Col<int> col_idx_a(nnz, arma::fill::zeros);
      arma::Col<int> col_idx_b(nnz, arma::fill::zeros);

      h5utl->read_sparse_coo_opdm(val_a, row_idx_a, col_idx_a, is_active, "D1A");
      h5utl->read_sparse_coo_opdm(val_b, row_idx_b, col_idx_b, is_active, "D1B");
      delete h5utl;

      for (size_t n = 0; n < nnz; n++) {
	 size_t i_a = row_idx_a(n);
	 size_t i_b = row_idx_b(n);
	 size_t j_a = col_idx_a(n);
	 size_t j_b = col_idx_b(n);
	 /* since we stored upper triangular part of the symmetric 1-RDM,
	    we need a factor of 2.0 for off-diagonal elements here */
	 double s_a = (i_a != j_a) ? 2.0 : 1.0;
	 double s_b = (i_b != j_b) ? 2.0 : 1.0;
	 double tmp_val_a = val_a(n);
	 double tmp_val_b = val_b(n);
         for(size_t p = 0; p < npts; p++) {
            rho_a_x(p) += s_a * ( phi_x(i_a, p) * phi(j_a, p) + phi(i_a, p) * phi_x(j_a, p) ) * tmp_val_a;
            rho_b_x(p) += s_b * ( phi_x(i_b, p) * phi(j_b, p) + phi(i_b, p) * phi_x(j_b, p) ) * tmp_val_b;
            rho_a_y(p) += s_a * ( phi_y(i_a, p) * phi(j_a, p) + phi(i_a, p) * phi_y(j_a, p) ) * tmp_val_a;
            rho_b_y(p) += s_b * ( phi_y(i_b, p) * phi(j_b, p) + phi(i_b, p) * phi_y(j_b, p) ) * tmp_val_b;
            rho_a_z(p) += s_a * ( phi_z(i_a, p) * phi(j_a, p) + phi(i_a, p) * phi_z(j_a, p) ) * tmp_val_a;
            rho_b_z(p) += s_b * ( phi_z(i_b, p) * phi(j_b, p) + phi(i_b, p) * phi_z(j_b, p) ) * tmp_val_b;
	 }
      }
      set_rhoa_x(rho_a_x);
      set_rhob_x(rho_b_x);
      set_rhoa_y(rho_a_y);
      set_rhob_y(rho_b_y);
      set_rhoa_z(rho_a_z);
      set_rhob_z(rho_b_z);
      // rho_a_x.print();

      for (int p = 0; p < npts; p++) {
         sigma_aa(p) = ( rho_a_x(p) * rho_a_x(p) ) +  ( rho_a_y(p) * rho_a_y(p) ) + ( rho_a_z(p) * rho_a_z(p) );
         sigma_bb(p) = ( rho_b_x(p) * rho_b_x(p) ) +  ( rho_b_y(p) * rho_b_y(p) ) + ( rho_b_z(p) * rho_b_z(p) );
         sigma_ab(p) = ( rho_a_x(p) * rho_b_x(p) ) +  ( rho_a_y(p) * rho_b_y(p) ) + ( rho_a_z(p) * rho_b_z(p) );
      }
      set_sigma_aa(sigma_aa);
      set_sigma_ab(sigma_ab);
      set_sigma_bb(sigma_bb);
   } 

#if 0
   void MCPDFT::build_density_gradients() {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());
      int p{0}, nu{0}, sigma{0};
      double tmp_d1a{0.0}, tmp_d1b{0.0};
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
                           private(p, nu, sigma, \
			           tmp_d1a, tmp_d1b)
      {
         for (int sigma = 0; sigma < nbfs; sigma++) {
             for (int nu = 0; nu < nbfs; nu++) {
	        tmp_d1a = D1a(sigma, nu);
	        tmp_d1b = D1b(sigma, nu);
	        if (fabs(tmp_d1a) < tol && fabs(tmp_d1b) < tol )
	           continue;
                #pragma omp for schedule(static) 
                   for (int p = 0; p < npts; p++) {
                      rho_a_x(p) += ( phi_x(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_x(nu, p) ) * tmp_d1a;
                      rho_b_x(p) += ( phi_x(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_x(nu, p) ) * tmp_d1b;
                      rho_a_y(p) += ( phi_y(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_y(nu, p) ) * tmp_d1a;
                      rho_b_y(p) += ( phi_y(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_y(nu, p) ) * tmp_d1b;
                      rho_a_z(p) += ( phi_z(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_z(nu, p) ) * tmp_d1a;
                      rho_b_z(p) += ( phi_z(sigma, p) * phi(nu, p) + phi(sigma, p) * phi_z(nu, p) ) * tmp_d1b;
                      sigma_aa(p) += ( rho_a_x(p) * rho_a_x(p) ) +  ( rho_a_y(p) * rho_a_y(p) ) + ( rho_a_z(p) * rho_a_z(p) );
                      sigma_bb(p) += ( rho_b_x(p) * rho_b_x(p) ) +  ( rho_b_y(p) * rho_b_y(p) ) + ( rho_b_z(p) * rho_b_z(p) );
                      sigma_ab(p) += ( rho_a_x(p) * rho_b_x(p) ) +  ( rho_a_y(p) * rho_b_y(p) ) + ( rho_a_z(p) * rho_b_z(p) );
                   }
            }
         }
      }
      set_rhoa_x(rho_a_x);
      set_rhob_x(rho_b_x);
      set_rhoa_y(rho_a_y);
      set_rhob_y(rho_b_y);
      set_rhoa_z(rho_a_z);
      set_rhob_z(rho_b_z);

      for (int p = 0; p < npts; p++) {
         sigma_aa(p) = ( rho_a_x(p) * rho_a_x(p) ) +  ( rho_a_y(p) * rho_a_y(p) ) + ( rho_a_z(p) * rho_a_z(p) );
         sigma_bb(p) = ( rho_b_x(p) * rho_b_x(p) ) +  ( rho_b_y(p) * rho_b_y(p) ) + ( rho_b_z(p) * rho_b_z(p) );
         sigma_ab(p) = ( rho_a_x(p) * rho_b_x(p) ) +  ( rho_a_y(p) * rho_b_y(p) ) + ( rho_a_z(p) * rho_b_z(p) );
      }
      set_sigma_aa(sigma_aa);
      set_sigma_ab(sigma_ab);
      set_sigma_bb(sigma_bb);
   }
#endif

   void MCPDFT::build_pi(const arma::mat &D2ab) {
      const bool is_sparse = MCPDFT::is_sparse();

      if(!is_sparse) { /* using dens RDMs*/
         build_ontop_pair_density(D2ab);
         if ( is_gga_ ) {
            build_ontop_pair_density_gradients(D2ab);
         }
      }else{ /* using sparse RDMs*/
         build_sparse_ontop_pair_density();
         if ( is_gga_ ) {
            build_sparse_ontop_pair_density_gradients();
         }
      }
   }

   void MCPDFT::build_ontop_pair_density(const arma::mat &D2ab) {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
      size_t npts = get_npts();
      arma::vec temp(npts, arma::fill::zeros);
      arma::mat phi(get_phi());
      int p{0},
	  mu{0}, nu{0},
	  lambda{0}, sigma{0};
      double tmp_d2ab{0.0};
      // pi(r,r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
      #pragma omp parallel default(shared) \
                           private(mu,nu,lambda,sigma, tmp_d2ab)
      {
         for (int mu = 0; mu < nbfs; mu++) {
            for (int nu = 0; nu < nbfs; nu++) {
               for (int lambda = 0; lambda < nbfs; lambda++) {
                  for (int sigma = 0; sigma < nbfs; sigma++) {
                     tmp_d2ab = D2ab(lambda*nbfs+sigma, mu*nbfs+nu);
            	     if (fabs(tmp_d2ab) < tol)
                        continue;
                        #pragma omp for schedule(static) 
                           for (int p = 0; p < npts; p++) {
                              temp(p) += phi(mu, p) * phi(nu, p) * phi(lambda, p) * phi(sigma, p) * tmp_d2ab;
            	           }
                  }
	       }
            }
         }
      } /* end of omp parallel for loop */
      set_pi(temp);
//      temp.print();
   }

   void MCPDFT::build_sparse_ontop_pair_density() {
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::vec temp(npts, arma::fill::zeros);
      bool is_active = MCPDFT::is_active();
      size_t nnz{0};
      HDF5Utility* h5utl = new HDF5Utility();

      h5utl->read_nnz(nnz, is_active, "D2AB");
      //std::printf("nnz = %6.ld\n",nnz);
      
      arma::vec val_ab(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim1(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim2(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim3(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim4(nnz, arma::fill::zeros);

      h5utl->read_sparse_coo_tpdm(val_ab,
		                  idx_dim1,
		                  idx_dim2,
		                  idx_dim3,
		                  idx_dim4,
				  is_active,
				  "D2AB");
      //val_ab.print();
      delete h5utl;

      for (size_t n = 0; n < nnz; n++) {
	 size_t idx1 = idx_dim1(n);
	 size_t idx2 = idx_dim2(n);
	 size_t idx3 = idx_dim3(n);
	 size_t idx4 = idx_dim4(n);
         /* since we stored upper triangular part of the symmetric 2-RDM,
         we need a factor of 2.0 for off-diagonal elements here */
	 double s_ab = (idx1 != idx2 || idx3 != idx4) ? 2.0 : 1.0 ;
	 double tmp_val_ab = val_ab(n);
         for(size_t p = 0; p < npts; p++) {
            temp(p) += s_ab * phi(idx1, p) * phi(idx2, p) * phi(idx3, p) * phi(idx4, p) * tmp_val_ab;
	 }
      }
      set_pi(temp);
//      temp.print();
   }
#if 0   
   void MCPDFT::build_ontop_pair_density(const arma::mat &D2ab) {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
      size_t npts = get_npts();
      arma::vec temp(npts, arma::fill::zeros);
      arma::mat phi(get_phi());
      int p{0},
	  mu{0}, nu{0},
	  lambda{0}, sigma{0};
      double tmp_d2ab{0.0};
      #pragma omp parallel default(shared) \
                     private(p, mu, nu, lambda, sigma)
      {
         // pi(r,r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
         for (int mu = 0; mu < nbfs; mu++) {
             for (int nu = 0; nu < nbfs; nu++) {
                 for (int lambda = 0; lambda < nbfs; lambda++) {
                     for (int sigma = 0; sigma < nbfs; sigma++) {
			/* careful since HDF5_READ transposes the matrix
			   it should not matter as long as it is symmetric */
                        tmp_d2ab = D2ab(lambda*nbfs+sigma, mu*nbfs+nu);
 			if (fabs(tmp_d2ab) < tol)
                           continue;
                        #pragma omp for schedule(static)
                           for (int p = 0; p < npts; p++) {
                              temp(p) += phi(mu, p) * phi(nu, p) * phi(lambda, p) * phi(sigma, p) * tmp_d2ab;
            	           }
                     }
                 }
             }
         } /* end of omp parallel for loop */
      }
      set_pi(temp);
//      temp.print();
//      double val = arma::sum(temp);
//      std::printf("sum of all elements: %-10.12lf\n", val);
   }
#endif

   void MCPDFT::build_ontop_pair_density_gradients(const arma::mat &D2ab) {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
      size_t npts = get_npts();
      int p{0},
	  mu{0}, nu{0},
	  lambda{0}, sigma{0};
      double tmp_d2ab{0.0};
      arma::mat phi(get_phi());
      arma::mat phi_x(get_phi_x());
      arma::mat phi_y(get_phi_y());
      arma::mat phi_z(get_phi_z());
      arma::vec pi_x(npts, arma::fill::zeros);
      arma::vec pi_y(npts, arma::fill::zeros);
      arma::vec pi_z(npts, arma::fill::zeros);
      #pragma omp parallel default(shared) \
                           private(p, mu, nu, lambda, sigma, tmp_d2ab)
      {
         // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
         for (int mu = 0; mu < nbfs; mu++) {
             for (int nu = 0; nu < nbfs; nu++) {
                 for (int lambda = 0; lambda < nbfs; lambda++) {
                     for (int sigma = 0; sigma < nbfs; sigma++) {
                        tmp_d2ab = D2ab(lambda*nbfs+sigma, mu*nbfs+nu);
                        if (fabs(tmp_d2ab) < tol)
                           continue;
                        #pragma omp for schedule(static)
                           for (int p = 0; p < npts; p++) {
                              pi_x(p) += ( phi_x(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi_x(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi_x(sigma, p)  * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_x(nu, p) ) * tmp_d2ab;

                              pi_y(p) += ( phi_y(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi_y(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi_y(sigma, p)  * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_y(nu, p) ) * tmp_d2ab;

                              pi_z(p) += ( phi_z(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi_z(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi_z(sigma, p)  * phi(nu, p) +
                                           phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_z(nu, p) ) * tmp_d2ab;
	 	          }
                     }
                 }
             }
         }
      }
      set_pi_x(pi_x);
      set_pi_z(pi_y);
      set_pi_z(pi_z);
   }

   void MCPDFT::build_sparse_ontop_pair_density_gradients() {
      size_t npts = get_npts();
      arma::mat phi(get_phi());
      arma::mat phi_x(get_phi_x());
      arma::mat phi_y(get_phi_y());
      arma::mat phi_z(get_phi_z());
      arma::vec pi_x(npts, arma::fill::zeros);
      arma::vec pi_y(npts, arma::fill::zeros);
      arma::vec pi_z(npts, arma::fill::zeros);
      arma::vec temp(npts, arma::fill::zeros);
      bool is_active = MCPDFT::is_active();
      size_t nnz{0};
      HDF5Utility* h5utl = new HDF5Utility();

      h5utl->read_nnz(nnz, is_active, "D2AB");
      //std::printf("nnz = %6.ld\n",nnz);
      
      arma::vec val_ab(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim1(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim2(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim3(nnz, arma::fill::zeros);
      arma::Col<int> idx_dim4(nnz, arma::fill::zeros);

      h5utl->read_sparse_coo_tpdm(val_ab,
		                  idx_dim1,
		                  idx_dim2,
		                  idx_dim3,
		                  idx_dim4,
				  is_active,
				  "D2AB");
      //val_ab.print();
      delete h5utl;

      for (size_t n = 0; n < nnz; n++) {
	 size_t idx1 = idx_dim1(n);
	 size_t idx2 = idx_dim2(n);
	 size_t idx3 = idx_dim3(n);
	 size_t idx4 = idx_dim4(n);
         /* since we stored upper triangular part of the symmetric 2-RDM,
         we need a factor of 2.0 for off-diagonal elements here */
	 double s_ab = (idx1 != idx2 || idx3 != idx4) ? 2.0 : 1.0 ;
	 double tmp_val_ab = val_ab(n);
         for(size_t p = 0; p < npts; p++) {
            pi_x(p) += s_ab * ( phi_x(idx1, p) * phi(idx2, p)   * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi_x(idx2, p) * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi_x(idx3, p)  * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi(idx3, p)    * phi_x(idx4, p) ) * tmp_val_ab;

            pi_y(p) += s_ab * ( phi_y(idx1, p) * phi(idx2, p)   * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi_y(idx2, p) * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi_y(idx3, p)  * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi(idx3, p)    * phi_y(idx4, p) ) * tmp_val_ab;

            pi_z(p) += s_ab * ( phi_z(idx1, p) * phi(idx2, p)   * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi_z(idx2, p) * phi(idx3, p)    * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi_z(idx3, p)  * phi(idx4, p) +
                                phi(idx1, p)   * phi(idx2, p)   * phi(idx3, p)    * phi_z(idx4, p) ) * tmp_val_ab;
	 }
      }
      set_pi_x(pi_x);
      set_pi_z(pi_y);
      set_pi_z(pi_z);
   }

#if 0
   void MCPDFT::build_ontop_pair_density_gradients(const arma::mat &D2ab) {
      double tol = 1.0e-20;
      int nbfs = get_nmo();
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
      double tmp_d2ab{0.0};
      #pragma omp parallel default(shared) \
                           private(p, mu, nu, lambda, sigma)
      {
         // pi(r,r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
         for (int mu = 0; mu < nbfs; mu++) {
             for (int nu = 0; nu < nbfs; nu++) {
                 for (int lambda = 0; lambda < nbfs; lambda++) {
                     for (int sigma = 0; sigma < nbfs; sigma++) {
			/* careful since HDF5_READ transposes the matrix
			   it should not matter as long as it is symmetric */
                        tmp_d2ab = D2ab(lambda*nbfs+sigma, mu*nbfs+nu);
 			if (fabs(tmp_d2ab) < tol)
                           continue;
                        #pragma omp for schedule(static)
                           for (int p = 0; p < npts; p++) {
                               pi_x(p) += ( phi_x(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi_x(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi_x(sigma, p)  * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_x(nu, p) ) * tmp_d2ab;

                               pi_y(p) += ( phi_y(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi_y(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi_y(sigma, p)  * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_y(nu, p) ) * tmp_d2ab;

                               pi_z(p) += ( phi_z(mu, p) * phi(lambda, p)   * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi_z(lambda, p) * phi(sigma, p)    * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi_z(sigma, p)  * phi(nu, p) +
                                            phi(mu, p)   * phi(lambda, p)   * phi(sigma, p)    * phi_z(nu, p) ) * tmp_d2ab;
            	           }
                     }
                 }
             }
         } /* end of omp parallel for loop */
      }
      set_pi_x(pi_x);
      set_pi_z(pi_y);
      set_pi_z(pi_z);
      
//      pi_x.print();
//      double val = arma::sum(pi_x);
//      std::printf("sum of all elements: %-10.12lf\n", val);
   }
#endif

   void MCPDFT::build_R() {
        double tol = 1.0e-20;
        size_t npts = get_npts();
        arma::vec pi(get_pi());
        arma::vec rho(get_rho());
        arma::vec temp(npts, arma::fill::zeros);
	double tmp_rho{0.0};
        #pragma parallel for schedule(static) \
	                     private(p) \
                             shared(npts, pi, rho, temp)
           for (int p = 0; p < npts; p++) {
              tmp_rho = rho(p);
              if ( tmp_rho < tol)
                 continue;
              temp(p) = 4.0 * pi(p) / ( tmp_rho * tmp_rho );
           }
        set_R(temp);
   }

#if 0
   void MCPDFT::build_R() {
        double tol = 1.0e-20;
        size_t npts = get_npts();
        arma::vec pi(get_pi());
        arma::vec rho(get_rho());
        arma::vec temp(npts, arma::fill::zeros);
	double tmp_rho{0.0};
        #pragma parallel for schedule(static) \
	                     private(p) \
                             shared(npts, pi, rho, temp)
           for (int p = 0; p < npts; p++) {
	       tmp_rho = rho(p);
	       if ( tmp_rho < tol)
	          continue;
               temp(p) = 4.0 * pi(p) / ( tmp_rho * tmp_rho );
           }
        set_R(temp);
   }
#endif

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
      arma::vec tr_rhoa(npts, arma::fill::zeros);
      arma::vec tr_rhob(npts, arma::fill::zeros);
      arma::vec tr_rho(npts, arma::fill::zeros);
      double rho = 0.0;
      double pi = 0.0;
      double zeta = 0.0;
      double R = 0.0;
      #pragma parallel for schedule(static) \
                           default(shared) \
                           private(p, zeta, pi, rho, R) \
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
         }
      set_tr_rhoa(tr_rhoa);
      set_tr_rhob(tr_rhob);
      
      tr_rho = tr_rhoa + tr_rhob;
      set_tr_rho(tr_rho);

      double dum_a = arma::dot(tr_rhoa, W);
      double dum_b = arma::dot(tr_rhob, W);
      double dum_tot = arma::dot(tr_rho, W);
      printf("\n");
      printf("  Integrated translated alpha density = %20.12lf\n",dum_a);
      printf("  Integrated translated beta density  = %20.12lf\n",dum_b);
      printf("  Integrated translated total density = %20.12lf\n",dum_tot);
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
     #pragma parallel for schedule(static) \
                          default(shared) \
                          private(p, zeta, pi, rho, R,\
		                  rho_x, rho_y, rho_z)
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

   void MCPDFT::fully_translate() {
     fully_translate_density();
     if ( is_gga_ ) {
        fully_translate_density_gradients();
     }
   }

   void MCPDFT::fully_translate_density(){
        double tol = 1.0e-20;
        double const R0 = 0.9;
        double const R1 = 1.15;
        double const A = -475.60656009;
        double const B = -379.47331922;
        double const C = -85.38149682;
        size_t npts = get_npts();
        arma::vec rho_vec(get_rho());
        arma::vec pi_vec(get_pi());
        arma::vec R_vec(get_R());
        arma::vec W(get_w());
        arma::vec tr_rhoa(npts, arma::fill::zeros);
        arma::vec tr_rhob(npts, arma::fill::zeros);
        arma::vec tr_rho(npts, arma::fill::zeros);
        double rho = 0.0;
        double pi = 0.0;
        double zeta = 0.0;
        double R = 0.0;
        double DelR = 0.0;
        #pragma parallel for schedule(static) \
                             default(shared) \
                             private(p, zeta, pi, rho) \
                             reduction(+:dum_a, dum_b, dum_tot)
        for (int p = 0; p < npts; p++) {
            zeta = 0.0;
            R = 0.0;
            rho = rho_vec(p);
            pi = pi_vec(p);
            DelR = R_vec(p) - R1;
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
        }
        set_tr_rhoa(tr_rhoa);
        set_tr_rhob(tr_rhob);

        tr_rho = tr_rhoa + tr_rhob;
        set_tr_rho(tr_rho);

        double temp_a = arma::dot(tr_rhoa, W);
        double temp_b = arma::dot(tr_rhob, W);
        double temp_tot = arma::dot(tr_rho, W);
        printf("\n");
        printf("      Integrated fully translated total density = %20.12lf\n",temp_tot);
        printf("      Integrated fully translated alpha density = %20.12lf\n",temp_a);
        printf("      Integrated fully translated beta density  = %20.12lf\n",temp_b);
        printf("\n");
   }

   void MCPDFT::fully_translate_density_gradients(){
        double tol = 1.0e-20;
        double const R0 = 0.9;
        double const R1 = 1.15;
        double const A = -475.60656009;
        double const B = -379.47331922;
        double const C = -85.38149682;
        size_t npts = get_npts();
        arma::vec rho_vec(get_rho());
        arma::vec pi_vec(get_pi());
        arma::vec R_vec(get_R());
        arma::vec W(get_w());
        arma::vec tr_rho_a_x(npts,  arma::fill::zeros);
        arma::vec tr_rho_b_x(npts,  arma::fill::zeros);
        arma::vec tr_rho_a_y(npts,  arma::fill::zeros);
        arma::vec tr_rho_b_y(npts,  arma::fill::zeros);
        arma::vec tr_rho_a_z(npts,  arma::fill::zeros);
        arma::vec tr_rho_b_z(npts,  arma::fill::zeros);
        arma::vec tr_sigma_aa(npts, arma::fill::zeros);
        arma::vec tr_sigma_ab(npts, arma::fill::zeros);
        arma::vec tr_sigma_bb(npts, arma::fill::zeros);
        double rho_x = 0.0;
        double rho_y = 0.0;
        double rho_z = 0.0;
        double rho = 0.0; 
        double pi = 0.0;
        double DelR = 0.0;
        double zeta = 0.0;
        double R = 0.0;
        #pragma parallel for schedule(static) \
                             default(shared) \
                             private(p, zeta, pi, rho, R, DelR,\
                                     rho_x, rho_y, rho_z)
           for (int p = 0; p < npts; p++) {
               rho_x = rho_a_x_(p) + rho_b_x_(p);
               rho_y = rho_a_y_(p) + rho_b_y_(p);
               rho_z = rho_a_z_(p) + rho_b_z_(p);
               rho = rho_vec(p);
               pi = pi_vec(p);
               DelR = R_vec(p) - R1;
               zeta = 0.0;
               R = 0.0;
               if ( !(rho < tol) && !(pi < 0.0) ) {
                   R = R_vec(p);
                   if ( ((1.0 - R) > tol) && ( R < R0 ) ) {
                      zeta = sqrt(1.0 - R);
                      tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0) + (R * rho_x) / (2.0*zeta) - pi_x_(p) / (rho*zeta);
                      tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0) - (R * rho_x) / (2.0*zeta) + pi_x_(p) / (rho*zeta);
                      tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0) + (R * rho_y) / (2.0*zeta) - pi_y_(p) / (rho*zeta);
                      tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0) - (R * rho_y) / (2.0*zeta) + pi_y_(p) / (rho*zeta);
                      tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0) + (R * rho_z) / (2.0*zeta) - pi_z_(p) / (rho*zeta);
                      tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0) - (R * rho_z) / (2.0*zeta) + pi_z_(p) / (rho*zeta);
                   }else if( !(R < R0) && !(R > R1) ) {
                           zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

                           tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0)
                                          + (A * pow(DelR, 4.0)) * ( (10.0 * pi_x_(p) / rho) - (5.0 * R * rho_x) )
                                          + (B * pow(DelR, 3.0)) * ( (8.0  * pi_x_(p) / rho) - (4.0 * R * rho_x) )
                                          + (C * pow(DelR, 2.0)) * ( (6.0  * pi_x_(p) / rho) - (3.0 * R * rho_x) );

                           tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0)
                                          + (A * pow(DelR, 4.0)) * (-(10.0 * pi_x_(p) / rho) + (5.0 * R * rho_x) )
                                          + (B * pow(DelR, 3.0)) * (-(8.0  * pi_x_(p) / rho) + (4.0 * R * rho_x) )
                                          + (C * pow(DelR, 2.0)) * (-(6.0  * pi_x_(p) / rho) + (3.0 * R * rho_x) );

                           tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0)
                                          + (A * pow(DelR, 4.0)) * ( (10.0 * pi_y_(p) / rho) - (5.0 * R * rho_y) )
                                          + (B * pow(DelR, 3.0)) * ( (8.0  * pi_y_(p) / rho) - (4.0 * R * rho_y) )
                                          + (C * pow(DelR, 2.0)) * ( (6.0  * pi_y_(p) / rho) - (3.0 * R * rho_y) );

                           tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0)
                                          + (A * pow(DelR, 4.0)) * (-(10.0 * pi_y_(p) / rho) + (5.0 * R * rho_y) )
                                          + (B * pow(DelR, 3.0)) * (-(8.0  * pi_y_(p) / rho) + (4.0 * R * rho_y) )
                                          + (C * pow(DelR, 2.0)) * (-(6.0  * pi_y_(p) / rho) + (3.0 * R * rho_y) );

                           tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0)
                                          + (A * pow(DelR, 4.0)) * ( (10.0 * pi_z_(p) / rho) - (5.0 * R * rho_z) )
                                          + (B * pow(DelR, 3.0)) * ( (8.0  * pi_z_(p) / rho) - (4.0 * R * rho_z) )
                                          + (C * pow(DelR, 2.0)) * ( (6.0  * pi_z_(p) / rho) - (3.0 * R * rho_z) );

                           tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0)
                                           + (A * pow(DelR, 4.0)) * (-(10.0 * pi_z_(p) / rho) + (5.0 * R * rho_z) )
                                           + (B * pow(DelR, 3.0)) * (-(8.0  * pi_z_(p) / rho) + (4.0 * R * rho_z) )
                                           + (C * pow(DelR, 2.0)) * (-(6.0  * pi_z_(p) / rho) + (3.0 * R * rho_z) );
                   }else if( R > R1 ) {
                           zeta = 0.0;
                           tr_rho_a_x(p) = (1.0 + zeta) * (rho_x/2.0);
                           tr_rho_b_x(p) = (1.0 - zeta) * (rho_x/2.0);
                           tr_rho_a_y(p) = (1.0 + zeta) * (rho_y/2.0);
                           tr_rho_b_y(p) = (1.0 - zeta) * (rho_y/2.0);
                           tr_rho_a_z(p) = (1.0 + zeta) * (rho_z/2.0);
                           tr_rho_b_z(p) = (1.0 - zeta) * (rho_z/2.0);
                   }
               }else{
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

   double MCPDFT::Hartree_energy(const arma::mat &D1a,
		                 const arma::mat &D1b,
		                 const arma::mat &Ja,
		                 const arma::mat &Jb) {
      double tmp{0.0};
      arma::mat dum(D1a);

      dum = arma::cdot(D1a,Ja);
      tmp  = arma::trace(dum);

      dum = arma::cdot(D1a,Jb);
      tmp += arma::trace(dum);

      dum = arma::cdot(D1b,Ja);
      tmp += arma::trace(dum);

      dum = arma::cdot(D1b,Jb);
      tmp += arma::trace(dum);

      tmp *= 0.5;
      return tmp;
   }
   
   double MCPDFT::core_energy(const arma::mat &D1,
		              const arma::mat &Hcore) {
      double tmp{0.0};
      arma::mat dum(D1);
      dum = arma::cdot(D1,Hcore);
      tmp = arma::trace(dum);
      return tmp;
   }
 
   bool MCPDFT::is_ao() const { return is_ao_; }
   bool MCPDFT::is_gga() const { return is_gga_; }
   bool MCPDFT::is_sparse() const { return is_sparse_; }
   bool MCPDFT::is_active() const { return is_active_; }
   size_t MCPDFT::get_nnz() const { return nnz_; }
   size_t MCPDFT::get_npts() const { return npts_; }
   size_t MCPDFT::get_nbfs() const { return nbfs_; }
   size_t MCPDFT::get_nao() const { return nao_; }
   size_t MCPDFT::get_nmo() const { return nmo_; }
   size_t MCPDFT::get_nactele() const { return nactele_; }
   size_t MCPDFT::get_nactorb() const { return nactorb_; }
   size_t MCPDFT::get_ncore() const { return ncore_; }
   size_t MCPDFT::get_nfrz() const { return nfrz_; }
   size_t MCPDFT::get_nocc() const { return nocc_; }
   size_t MCPDFT::get_nvir() const { return nvir_; }
   arma::vec MCPDFT::get_w() const { return w_; }
   arma::vec MCPDFT::get_x() const { return x_; }
   arma::vec MCPDFT::get_y() const { return y_; }
   arma::vec MCPDFT::get_z() const { return z_; }
   arma::mat MCPDFT::get_phi() const { return phi_; }
   arma::mat MCPDFT::get_phi_x() const { return phi_x_; }
   arma::mat MCPDFT::get_phi_y() const { return phi_y_; }
   arma::mat MCPDFT::get_phi_z() const { return phi_z_; }
   double MCPDFT::get_enuc() const { return enuc_; }
   double MCPDFT::get_eref() const { return eref_; }
   double MCPDFT::get_eclass()  const { return eclass_; }
   arma::mat MCPDFT::get_cmat() const { return cmat_; }
   arma::mat MCPDFT::get_hcore() const { return Hcore_; }
   arma::mat MCPDFT::get_ja() const { return ja_; }
   arma::mat MCPDFT::get_jb() const { return jb_; }
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

   void MCPDFT::set_nnz(const size_t nnz) { nnz_ = nnz; }
   void MCPDFT::set_npts(const size_t npts) { npts_ = npts; }
   void MCPDFT::set_nbfs(const size_t nbfs)    { nbfs_ = nbfs; }
   void MCPDFT::set_nao(const size_t nao) { nao_ = nao; }
   void MCPDFT::set_nmo(const size_t nmo) { nmo_ = nmo; }
   void MCPDFT::set_nactele(const size_t nactele) { nactele_ = nactele; }
   void MCPDFT::set_nactorb(const size_t nactorb) { nactorb_ = nactorb; }
   void MCPDFT::set_ncore(const size_t ncore) { ncore_ = ncore; }
   void MCPDFT::set_nfrz(const size_t nfrz) { nfrz_ = nfrz; }
   void MCPDFT::set_nocc(const size_t nocc) { nocc_ = nocc; }
   void MCPDFT::set_nvir(const size_t nvir) { nvir_ = nvir; }
   void MCPDFT::set_w(const arma::vec &w) { w_ = w; }
   void MCPDFT::set_x(const arma::vec &x) { x_ = x; }
   void MCPDFT::set_y(const arma::vec &y) { y_ = y; }
   void MCPDFT::set_z(const arma::vec &z) { z_ = z; }
   void MCPDFT::set_phi(const arma::mat &phi) { phi_ = phi; }
   void MCPDFT::set_phi_x(const arma::mat &phi_x) { phi_x_ = phi_x; }
   void MCPDFT::set_phi_y(const arma::mat &phi_y) { phi_y_ = phi_y; }
   void MCPDFT::set_phi_z(const arma::mat &phi_z) { phi_z_ = phi_z; }
   void MCPDFT::set_enuc(const double enuc) { enuc_ = enuc; }
   void MCPDFT::set_eref(const double eref) { eref_ = eref; }
   void MCPDFT::set_eclass(const double eclass) { eclass_ = eclass; }
   void MCPDFT::set_cmat(const arma::mat &cmat) { cmat_ = cmat; }
   void MCPDFT::set_hcore(const arma::mat &Hcore) { Hcore_ = Hcore; }
   void MCPDFT::set_ja(const arma::mat &ja) { ja_ = ja; }
   void MCPDFT::set_jb(const arma::mat &jb) { jb_ = jb; }
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

      printf("   # M. Mostafanejad and A. E. DePrince III\n");
      printf("      J. Chem. Theory Comput. 15, 290-302 (2019).\n\n");
      printf("   # M. Mostafanejad, M. D. Liebenthal, and A. E. DePrince III\n");
      printf("      J. Chem. Theory Comput. 16, 2274-2283 (2020).\n\n");
   }
}
