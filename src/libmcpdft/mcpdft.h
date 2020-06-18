#ifndef MCPDFT_H
#define MCPDFT_H

#include <armadillo>

namespace mcpdft{

class MCPDFT {
   public:
      /// constructor
      MCPDFT();
      MCPDFT(std::string test_case);
      /// destructor
      ~MCPDFT();

      /// initialize the class member variables using UTF-format input files
      void common_init(std::string test_case);

      /// initialize the class member variables using HDF5 input files 
      void common_init();

      /// build the 1-particle density matrices (OPDMs)
      void build_opdm();

      /// build the alpha-beta block of 2-particle density matrix (TPDM)
      void build_tpdm();

      /// build density functions and their gradients on the grid points
      void build_rho();

      /// build spin and total density functions rhoa(r), rhob(r) and rho(r)
      void build_density_functions();

      /// build density gradients rhoa_x(r), rhoa_y(r) ... on the grid points
      void build_density_gradients();

      /// build on-top pair-density pi(r,r) and its gradinets pi_x(r,r) ... on the grids
      void build_pi(const arma::mat &D2ab);

      /// build on-top pair-density pi(r,r) on the grids
      void build_ontop_pair_density(const arma::mat &D2ab);

      /// build on-top pair-density gradinets pi_x(r,r), pi_y(r,r) ... on the grids
      void build_ontop_pair_density_gradients(const arma::mat &D2ab);

      /// build R(r) factor
      void build_R();

      /// translate the one-electron densities and their gradients
      void translate();

      /// translate the one-electron density functions
      void translate_density();

      /// translate the density gradients
      void translate_density_gradients();

      /// fully translate the one-electron densities 
      void fully_translate();

      /// fully translate the one-electron density functions
      void fully_translate_density();

      /// fully translate the density gradients
      void fully_translate_density_gradients();

      //=============== accessors ===============//
      bool is_gga() const;
      size_t get_npts() const;
      size_t get_nbfs() const;
      size_t get_nao() const;
      size_t get_nmo() const;
      size_t get_nactele() const;
      size_t get_nactorb() const;
      size_t get_ncore() const;
      size_t get_nfrz() const; 
      size_t get_nocc() const;
      size_t get_nvir() const;
      arma::vec get_w() const;
      arma::vec get_x() const;
      arma::vec get_y() const;
      arma::vec get_z() const;
      arma::mat get_phi() const;
      arma::mat get_phi_x() const;
      arma::mat get_phi_y() const;
      arma::mat get_phi_z() const;
      double get_eref() const;
      double get_eclass() const;
      arma::mat get_cmat() const;
      arma::mat get_hcore() const;
      arma::mat get_ja() const;
      arma::mat get_jb() const;
      arma::mat get_D1a() const;
      arma::mat get_D1b() const;
      arma::mat get_D2ab() const;
      arma::vec get_rhoa() const;
      arma::vec get_rhoa_x() const;
      arma::vec get_rhoa_y() const;
      arma::vec get_rhoa_z() const;
      arma::vec get_rhob() const;
      arma::vec get_rhob_x() const;
      arma::vec get_rhob_y() const;
      arma::vec get_rhob_z() const;
      arma::vec get_rho() const;
      arma::vec get_tr_rhoa() const;
      arma::vec get_tr_rhoa_x() const;
      arma::vec get_tr_rhoa_y() const;
      arma::vec get_tr_rhoa_z() const;
      arma::vec get_tr_rhob() const;
      arma::vec get_tr_rhob_x() const;
      arma::vec get_tr_rhob_y() const;
      arma::vec get_tr_rhob_z() const;
      arma::vec get_tr_rho() const;
      arma::vec get_pi() const;
      arma::vec get_R() const;
      arma::vec get_sigma_aa() const;
      arma::vec get_sigma_ab() const;
      arma::vec get_sigma_bb() const;
      arma::vec get_tr_sigma_aa() const;
      arma::vec get_tr_sigma_ab() const;
      arma::vec get_tr_sigma_bb() const;

      void set_npts(const size_t npts);
      void set_nbfs(const size_t nbfs);
      void set_nao(const size_t nao);
      void set_nmo(const size_t nmo);
      void set_nactele(const size_t nactele);
      void set_nactorb(const size_t nactorb);
      void set_ncore(const size_t ncore);
      void set_nfrz(const size_t nfrz);
      void set_nocc(const size_t nocc);
      void set_nvir(const size_t nvir);
      void set_w(const arma::vec &w);
      void set_x(const arma::vec &x);
      void set_y(const arma::vec &y);
      void set_z(const arma::vec &z);
      void set_phi(const arma::mat &phi);
      void set_phi_x(const arma::mat &phi_x);
      void set_phi_y(const arma::mat &phi_y);
      void set_phi_z(const arma::mat &phi_z);
      void set_eref(const double eref);
      void set_eclass(const double eclass);
      void set_cmat(const arma::mat &cmat);
      void set_hcore(const arma::mat &HCore);
      void set_ja(const arma::mat &ja);
      void set_jb(const arma::mat &jb);
      void set_D1a(const arma::mat &D1a);
      void set_D1b(const arma::mat &D1b);
      void set_D2ab(const arma::mat &D2ab);
      void set_rhoa(const arma::vec &rhoa);
      void set_rhoa_x(const arma::vec &rhoa_x);
      void set_rhoa_y(const arma::vec &rhoa_y);
      void set_rhoa_z(const arma::vec &rhoa_z);
      void set_rhob(const arma::vec &rhob);
      void set_rhob_x(const arma::vec &rhob_x);
      void set_rhob_y(const arma::vec &rhob_y);
      void set_rhob_z(const arma::vec &rhob_z);
      void set_rho(const arma::vec &rho);
      void set_tr_rhoa(const arma::vec &tr_rhoa);
      void set_tr_rhoa_x(const arma::vec &tr_rhoa_x);
      void set_tr_rhoa_y(const arma::vec &tr_rhoa_y);
      void set_tr_rhoa_z(const arma::vec &tr_rhoa_z);
      void set_tr_rhob(const arma::vec &tr_rhob);
      void set_tr_rhob_x(const arma::vec &tr_rhob_x);
      void set_tr_rhob_y(const arma::vec &tr_rhob_y);
      void set_tr_rhob_z(const arma::vec &tr_rhob_z);
      void set_tr_rho(const arma::vec &tr_rho);
      void set_pi(const arma::vec &pi);
      void set_pi_x(const arma::vec &pi_x);
      void set_pi_y(const arma::vec &pi_y);
      void set_pi_z(const arma::vec &pi_z);
      void set_R(const arma::vec &R);
      void set_sigma_aa(const arma::vec &sigma_aa);
      void set_sigma_ab(const arma::vec &sigma_ab);
      void set_sigma_bb(const arma::vec &sigma_bb);
      void set_tr_sigma_aa(const arma::vec &tr_sigma_aa);
      void set_tr_sigma_ab(const arma::vec &tr_sigma_ab);
      void set_tr_sigma_bb(const arma::vec &tr_sigma_bb);
      //==========================================// end of accessors
      
      void print_banner() const;

   protected:   
       //=========== utility functions ============//
       void read_grids_from_file(std::string test_case);
       void read_orbitals_from_file(std::string test_case);
       void read_gradients_from_file(std::string test_case);
       void read_energies_from_file(std::string test_case);
       void read_opdm_from_file(std::string test_case);
       void read_cmat_from_file(std::string test_case);
       //==========================================// end utility functions

   private:
       /// A boolean variable to show if a chose functional is GGA or not
       bool is_gga_;

       /// 1-electron reduced-density matrix of alpha spin
       arma::mat D1a_;
 
       /// 1-electron reduced-density matrix of beta spin
       arma::mat D1b_;
 
       /// 2-electron reduced-density matrix of alpha-beta spin
       arma::mat D2ab_;

      /// number of grid points
      size_t npts_;

      /// number of basis functions
      size_t nbfs_;
      
      /// number of atomic orbitals (AOs)
      size_t nao_;

      /// number of molecular orbitals (MOs)
      size_t nmo_;

      /// number of active electrons
      size_t nactele_;

      /// number of active orbitals
      size_t nactorb_;

      /// number of core orbitals
      size_t ncore_;

      /// number of frozen orbitals
      size_t nfrz_;

      /// number of occupied (inactive) orbotals
      size_t nocc_;

      /// number of virtual (secondary or unoccupied) orbitals
      size_t nvir_;

      /// vector of weights for quadrature grid points
      arma::vec w_;

      /// vector of x-coordinates for quadrature grid points
      arma::vec x_;

      /// vector of y-coordinates for quadrature grid points
      arma::vec y_;

      /// vector of z-coordinates for quadrature grid points
      arma::vec z_;

      /// matrix of grids([w],[x],[y],[z]) columns [x], [y], ...
      arma::mat grids_;

      /// matrix of orbital values calculated on the grid points phi(npts, nbfs)
      arma::mat phi_;

      /// matrix of orbital x-derivative values calculated on the grid points phi(npts, nbfs)
      arma::mat phi_x_;

      /// matrix of orbital y-derivative values calculated on the grid points phi(npts, nbfs)
      arma::mat phi_y_;

      /// matrix of orbital z-derivative values calculated on the grid points phi(npts, nbfs)
      arma::mat phi_z_;

      /// the reference electronic energy
      double eref_;

      /// the classical electronic energy 
      double eclass_;

      /// AO->MO transformation matrix C
      arma::mat cmat_;

      /// core Hamiltonian matrix
      arma::mat Hcore_;

      /// classical Coulomb (Hartree) alpha elerctron-electron interaction matrix J
      arma::mat ja_;

      /// classical Coulomb (Hartree) beta elerctron-electron interaction matrix J
      arma::mat jb_;

      /// alpha density vector rho_a(r)
      arma::vec rho_a_;

      /// x-derivative of the alpha density vector rho_a(r)
      arma::vec rho_a_x_;

      /// y-derivative of the alpha density vector rho_a(r)
      arma::vec rho_a_y_;

      /// z-derivative of the alpha density vector rho_a(r)
      arma::vec rho_a_z_;

      /// beta density vector rho_b(r)
      arma::vec rho_b_;

      /// x-derivative of the beta density vector rho_b(r)
      arma::vec rho_b_x_;

      /// y-derivative of the beta density vector rho_b(r)
      arma::vec rho_b_y_;

      /// z-derivative of the beta density vector rho_b(r)
      arma::vec rho_b_z_;

      /// total density vector rho_(r)
      arma::vec rho_;

      /// (fully-)translated alpha density vector tr_rho_a(r)
      arma::vec tr_rho_a_;

      /// x-derivative of the (fully-)translated alpha density vector tr_rho_a(r)
      arma::vec tr_rho_a_x_;

      /// y-derivative of the (fully-)translated alpha density vector tr_rho_a(r)
      arma::vec tr_rho_a_y_;

      /// z-derivative of the (fully-)translated alpha density vector tr_rho_a(r)
      arma::vec tr_rho_a_z_;

      /// (fully-)translated beta density vector tr_rho_b(r)
      arma::vec tr_rho_b_;

      /// x-derivative of the (fully-)translated beta density vector tr_rho_b(r)
      arma::vec tr_rho_b_x_;

      /// y-derivative of the (fully-)translated beta density vector tr_rho_b(r)
      arma::vec tr_rho_b_y_;

      /// z-derivative of the (fully-)translated beta density vector tr_rho_b(r)
      arma::vec tr_rho_b_z_;

      /// (fully-)translated total density vector tr_rho_(r)
      arma::vec tr_rho_;

      /// on-top pair-density vector pi_(r)
      arma::vec pi_;

      /// x-derivative of the on-top pair-density vector pi_(r)
      arma::vec pi_x_;

      /// y-derivative of the on-top pair-density vector pi_(r)
      arma::vec pi_y_;

      /// z-derivative of the on-top pair-density vector pi_(r)
      arma::vec pi_z_;

      /// R(r) factor
      arma::vec R_;

      /// sigma alpha-alpha vector sigma_aa(r)
      arma::vec sigma_aa_;

      /// sigma alpha-beta vector sigma_ab(r)
      arma::vec sigma_ab_;

      /// sigma beta-beta vector sigma_bb(r)
      arma::vec sigma_bb_;

      /// (fully-)translated sigma alpha-alpha vector tr_sigma_aa(r)
      arma::vec tr_sigma_aa_;

      /// (fully-)translated sigma alpha-beta vector tr_sigma_ab(r)
      arma::vec tr_sigma_ab_;

      /// (fully-)translated sigma beta-beta vector tr_sigma_bb(r)
      arma::vec tr_sigma_bb_;
};

}
#endif // MCPDFT_H
