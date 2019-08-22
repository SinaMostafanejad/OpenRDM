#include <armadillo>
#include <iostream>
#include <fstream>

#include "mcpdft.h"

namespace mcpdft {

   MCPDFT::MCPDFT()  {};
   MCPDFT::~MCPDFT() {};

   void MCPDFT::common_init() {

       read_grids_from_file();
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
               file >> w(i) >> x(i) >> y(i) >> z(i) >> phi(i);
            }
            set_w(w);
            set_x(x);
            set_y(y);
            set_z(z);
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

   size_t MCPDFT::get_npts() const { return npts_; }
   arma::vec MCPDFT::get_w() const { return w_; }
   arma::vec MCPDFT::get_x() const { return x_; }
   arma::vec MCPDFT::get_y() const { return y_; }
   arma::vec MCPDFT::get_z() const { return z_; }
   arma::vec MCPDFT::get_phi() const { return phi_; }
   double MCPDFT::get_eref() const { return eref_; }
   double MCPDFT::get_eclass() const { return eclass_; }
   arma::mat MCPDFT::get_cmat() const { return cmat_; }

   void MCPDFT::set_npts(const size_t npts) { npts_ = npts; }
   void MCPDFT::set_w(const arma::vec &w) { w_ = w; }
   void MCPDFT::set_x(const arma::vec &x) { x_ = x; }
   void MCPDFT::set_y(const arma::vec &y) { y_ = y; }
   void MCPDFT::set_z(const arma::vec &z) { z_ = z; }
   void MCPDFT::set_phi(const arma::vec &phi) { phi_ = phi; }
   void MCPDFT::set_eref(const double eref) { eref_ = eref; }
   void MCPDFT::set_eclass(const double eclass) { eclass_ = eclass; }
   void MCPDFT::set_cmat(const arma::mat &cmat) { cmat_ = cmat; }

}
