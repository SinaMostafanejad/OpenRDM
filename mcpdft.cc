#include <armadillo>
#include <iostream>
#include <fstream>

#include "mcpdft.h"

namespace mcpdft {

   MCPDFT::MCPDFT()  { common_init(); }
   MCPDFT::~MCPDFT() {}

   void MCPDFT::common_init() {

       read_grids_from_file();
       read_orbitals_from_file();
       read_energies_from_file(); 
       read_cmat_from_file();
   }

   arma::vec MCPDFT::build_rho(arma::mat &phi,
                               arma::mat &D1) {
       
       D1 * phi * phi;

       return D1;
   
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
              for (int i = 0; i < nbfs; i++) { 
                  D1a(mu, nu) = ca(mu, i) * ca(nu, i);
                  D1b(mu, nu) = cb(mu, i) * cb(nu, i);
              }
          }
      }
      // D1a.print("D1a = ");
      // D1b.print("D1b = ");
      set_D1a(D1a);
      set_D1b(D1b);
   }

   void MCPDFT::build_tpdm() {

      // fetching the number of basis functions
      int nbfs;
      nbfs = get_nbfs();

      arma::mat D1a(get_D1a());
      arma::mat D1b(get_D1b());

      arma::mat D2ab(nbfs*nbfs, nbfs*nbfs, arma::fill::zeros);
      D2ab = arma::kron(D1a,D1b);
      // D2ab.print("D2ab = ");
      set_D2ab(D2ab);

   }


   size_t MCPDFT::get_npts() const { return npts_; }
   int    MCPDFT::get_nbfs() const { return nbfs_; }
   arma::vec MCPDFT::get_w() const { return w_; }
   arma::vec MCPDFT::get_x() const { return x_; }
   arma::vec MCPDFT::get_y() const { return y_; }
   arma::vec MCPDFT::get_z() const { return z_; }
   arma::mat MCPDFT::get_phi() const { return phi_; }
   double MCPDFT::get_eref() const { return eref_; }
   double MCPDFT::get_eclass()  const { return eclass_; }
   arma::mat MCPDFT::get_cmat() const { return cmat_; }
   arma::mat MCPDFT::get_D1a()  const { return D1a_ ; }
   arma::mat MCPDFT::get_D1b()  const { return D1b_ ; }
   arma::mat MCPDFT::get_D2ab()  const { return D2ab_ ; }

   void MCPDFT::set_npts(const size_t npts) { npts_ = npts; }
   void MCPDFT::set_nbfs(const int nbfs)    { nbfs_ = nbfs; }
   void MCPDFT::set_w(const arma::vec &w) { w_ = w; }
   void MCPDFT::set_x(const arma::vec &x) { x_ = x; }
   void MCPDFT::set_y(const arma::vec &y) { y_ = y; }
   void MCPDFT::set_z(const arma::vec &z) { z_ = z; }
   void MCPDFT::set_phi(const arma::mat &phi) { phi_ = phi; }
   void MCPDFT::set_eref(const double eref) { eref_ = eref; }
   void MCPDFT::set_eclass(const double eclass) { eclass_ = eclass; }
   void MCPDFT::set_cmat(const arma::mat &cmat) { cmat_ = cmat; }
   void MCPDFT::set_D1a(const arma::mat &D1a) { D1a_ = D1a; }
   void MCPDFT::set_D1b(const arma::mat &D1b) { D1b_ = D1b; }
   void MCPDFT::set_D2ab(const arma::mat &D2ab) { D2ab_ = D2ab; }

}
