#include "mcpdft.h"

namespace mcpdft {

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

   void MCPDFT::read_opdm_from_file() {

       std::ifstream file;

       // reading the AO->MO transformation matrix C from file
       file.open("./data/opdm.txt");

       if (!file)
          std::cout << "Error opening file.\n";
       else {
            int cols{ 0 };
            int rows{ 0 };
            file >> rows >> cols;
            arma::mat D1a(rows, cols, arma::fill::zeros);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    file >> D1a(i,j);

            set_D1a(D1a);
            set_D1b(D1a);
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

}
