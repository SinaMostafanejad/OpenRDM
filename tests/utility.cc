#include "mcpdft.h"

namespace mcpdft {

   void MCPDFT::read_grids_from_file(std::string test_case) {
       std::string grids_fname = std::string("./")
	                       + test_case
			       + std::string("/grids.txt");
       std::ifstream file;

       file.open(grids_fname);

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

   void MCPDFT::read_orbitals_from_file(std::string test_case) {
       std::string orbs_fname  = std::string("./")
	                       + test_case
			       + std::string("/orbitals.txt");
       std::ifstream file;

       file.open(orbs_fname);

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

   void MCPDFT::read_gradients_from_file(std::string test_case) {
       std::string grads_fname  = std::string("./")
	                        + test_case
			        + std::string("/gradients.txt");
       std::ifstream file;

       file.open(grads_fname);

       if (!file)
          std::cout << "Error opening file.\n";
       else {
            size_t npts = 0;
            int    nbfs = 0;
            file >> npts >> nbfs;
            set_npts(npts);
            set_nbfs(nbfs);

            arma::mat phi_x(npts, nbfs, arma::fill::zeros);
            arma::mat phi_y(npts, nbfs, arma::fill::zeros);
            arma::mat phi_z(npts, nbfs, arma::fill::zeros);

            for (int p = 0; p < npts; p++)
                for (int mu = 0; mu < nbfs; mu++)
                    file >> phi_x(p, mu) >> phi_y(p, mu) >> phi_z(p, mu);

            set_phi_x(phi_x);
            set_phi_y(phi_y);
            set_phi_z(phi_z);
       }
       file.close();
   }


   void MCPDFT::read_energies_from_file(std::string test_case) {
       std::string eref_fname  = std::string("./")
	                       + test_case
			       + std::string("/eref.txt");

       std::string eclass_fname  = std::string("./")
	                         + test_case
			         + std::string("/eclass.txt");
       std::ifstream file;

       // reading the reference energy from file
       file.open(eref_fname);

       if (!file)
          std::cout << "Error opening file.\n";
       else {
            double eref = 0.0;
            file >> eref;
            set_eref(eref);
       }
       file.close();

       // reading the classical energy from file
       file.open(eclass_fname);

       if (!file)
          std::cout << "Error opening file.\n";
       else {
            double eclass = 0.0;
            file >> eclass;
            set_eclass(eclass);
       }
       file.close();
   }

   void MCPDFT::read_opdm_from_file(std::string test_case) {
       std::string opdm_fname  = std::string("./")
	                       + test_case
			       + std::string("/opdm.txt");
       std::ifstream file;

       // reading the AO->MO transformation matrix C from file
       file.open(opdm_fname);

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

   void MCPDFT::read_cmat_from_file(std::string test_case) {
       std::string cmat_fname  = std::string("tests/")
	                       + test_case
			       + std::string("/cmat.txt");
       std::ifstream file;

       // reading the AO->MO transformation matrix C from file
       file.open(cmat_fname);

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
