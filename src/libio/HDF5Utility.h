#ifndef HDF5UTILITY_H
#define HDF5UTILITY_H

#include <armadillo>

namespace mcpdft {

class HDF5Utility {
   public:
      /// constructor
      HDF5Utility();

      /// destructor
      ~HDF5Utility();

      /// read the number of grid point numbers and basis functions (AOs, MOs, NO, etc.)
      void read_nbfs(size_t &nao,
		     size_t &nmo,
		     size_t &npts);

      /// read active space variables (nactele, nactorb, ncore, etc.)
      void read_active_space_vars(size_t &nactele,
		                  size_t &nactorb,
				  size_t &ncore, 
				  size_t &nfrz, 
				  size_t &nocc, 
				  size_t &nvir);

      /// read AO2MO transformation matrix C
      void read_AO2MO_Cmat(arma::mat &Cmat);

      /// read classical nuclear repulsion energy
      void read_enuc(double &enuc);

      /// read classical Coulomb (Hartree) electron-electron interaction matrices Ja and Jb
      void read_Hartree_Jmats(arma::mat &Ja,
		              arma::mat &Jb,
		              bool is_ao = true);

      /// read core Hamiltonian matrix Hcore
      void read_hcore(arma::mat &Hcore,
		      bool is_ao = true);

      /// read grids w, x, y, z
      void read_grids(arma::vec &W,
		      arma::vec &X,
                      arma::vec &Y,
                      arma::vec &Z);

      /// read AOs/MOs and their gratients
      void read_superphi(arma::mat &phi,
		         arma::mat &phi_x,
			 arma::mat &phi_y,
			 arma::mat &phi_z,
			 bool is_ao = true);

};

}
#endif // HDF5UTILITY_H
