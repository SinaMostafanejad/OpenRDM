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

   void MCPDFT::read_grids_from_file() {

       std::ifstream file;

       file.open("./data/grids.txt");
       
       if (!file)
          std::cout << "Error opening file.\n";
       else {
            size_t npts = 0;
            file >> npts;
            set_npts(npts);
            arma::vec w(npts);
            arma::vec x(npts); 
            arma::vec y(npts); 
            arma::vec z(npts); 
            arma::vec phi(npts); 
            for (int i = 0; i < npts; i++){
               file >> w(i) >> x(i) >> y(i) >> z(i);
            }
            set_w(w);
            set_x(x);
            set_y(y);
            set_z(z);
            set_phi(phi);
       }
       file.close();
   }

   void MCPDFT::read_orbitals_from_file() {

       std::ifstream file;

       file.open("./data/orbitals.txt");
       
       if (!file)
          std::cout << "Error opening file.\n";
       else {
            size_t npts = 0;
            int    nbfs = 0;
            file >> npts >> nbfs;
            set_npts(npts);
            set_nbfs(nbfs);

            arma::mat phi(npts, nbfs, arma::fill::zeros);

            for (int p = 0; p < npts; p++)
                for (int mu = 0; mu < nbfs; mu++)
                    file >> phi(p, mu);
            
            set_phi(phi);
       }
       file.close();
   }

   void MCPDFT::read_energies_from_file() {

       std::ifstream file;

       // reading the reference energy from file
       file.open("./data/eref.txt");
       
       if (!file)
          std::cout << "Error opening file.\n";
       else {
            double eref = 0.0;
            file >> eref;
            set_eref(eref);
       }
       file.close();

       // reading the classical energy from file
       file.open("./data/eclass.txt");
       
       if (!file)
          std::cout << "Error opening file.\n";
       else {
            double eclass = 0.0;
            file >> eclass;
            set_eclass(eclass);
       }
       file.close();
   }

   void MCPDFT::read_cmat_from_file() {

       std::ifstream file;

       // reading the AO->MO transformation matrix C from file
       file.open("./data/cmat.txt");
       
       if (!file)
          std::cout << "Error opening file.\n";
       else {
            int cols{ 0 };
            int rows{ 0 };
            file >> rows >> cols;
            arma::mat cmat(rows, cols, arma::fill::zeros);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    file >> cmat(i,j);

            set_cmat(cmat);
       }
       file.close();
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

}
