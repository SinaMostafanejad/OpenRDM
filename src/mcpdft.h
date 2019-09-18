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

      /// initialize the class member variables
      void common_init(std::string test_case);

      /// Computes the MCPDFT energy
      double mcpdft_energy(const arma::mat &D1a,
                           const arma::mat &D1b,
                           const arma::mat &D2ab);

      /// Build the 1-particle density matrices (OPDMs)
      void build_opdm();

      /// Build the alpha-beta block of 2-particle density matrix (TPDM)
      void build_tpdm();

      /// Build spin and total density functions rhoa(r), rhob(r) and rho(r)
      void build_rho();

      /// Build on-top pair-density pi(r,r)
      void build_pi(const arma::mat &D2ab);

      /// Build R(r) factor
      void build_R();

      /// Translate the one-electron densities 
      void translate();

      //=============== accessors ===============//
      size_t get_npts() const;
      int    get_nbfs() const;
      arma::vec get_w() const;
      arma::vec get_x() const;
      arma::vec get_y() const;
      arma::vec get_z() const;
      arma::mat get_phi() const;
      double get_eref() const;
      double get_eclass() const;
      arma::mat get_cmat() const;
      arma::mat get_D1a() const;
      arma::mat get_D1b() const;
      arma::mat get_D2ab() const;
      arma::vec get_rhoa() const;
      arma::vec get_rhob() const;
      arma::vec get_rho() const;
      arma::vec get_tr_rhoa() const;
      arma::vec get_tr_rhob() const;
      arma::vec get_tr_rho() const;
      arma::vec get_pi() const;
      arma::vec get_R() const;

      void set_npts(const size_t npts);
      void set_nbfs(const int nbfs);
      void set_w(const arma::vec &w);
      void set_x(const arma::vec &x);
      void set_y(const arma::vec &y);
      void set_z(const arma::vec &z);
      void set_phi(const arma::mat &phi);
      void set_eref(const double eref);
      void set_eclass(const double eclass);
      void set_cmat(const arma::mat &cmat);
      void set_D1a(const arma::mat &D1a);
      void set_D1b(const arma::mat &D1b);
      void set_D2ab(const arma::mat &D2ab);
      void set_rhoa(const arma::vec &rhoa);
      void set_rhob(const arma::vec &rhob);
      void set_rho(const arma::vec &rho);
      void set_tr_rhoa(const arma::vec &tr_rhoa);
      void set_tr_rhob(const arma::vec &tr_rhob);
      void set_tr_rho(const arma::vec &tr_rho);
      void set_pi(const arma::vec &pi);
      void set_R(const arma::vec &R);
      //==========================================// end of accessors

   protected:

      //=========== utility functions ============//
      void read_grids_from_file(std::string test_case);
      void read_orbitals_from_file(std::string test_case);
      void read_energies_from_file(std::string test_case);
      void read_opdm_from_file(std::string test_case);
      void read_cmat_from_file(std::string test_case);
      //==========================================// end of utility functions

   private:

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

      /// the reference electronic energy (read from file)
      double eref_;

      /// the classical electronic energy (read from file)
      double eclass_;

      /// AO->MO transformation matrix C (read from file)
      arma::mat cmat_;

      /// alpha density vector rho_a(r)
      arma::vec rho_a_;

      /// beta density vector rho_b(r)
      arma::vec rho_b_;

      /// total density vector rho_(r)
      arma::vec rho_;

      /// (fully-)translated alpha density vector tr_rho_a(r)
      arma::vec tr_rho_a_;

      /// (fully-)translated beta density vector tr_rho_b(r)
      arma::vec tr_rho_b_;

      /// (fully-)translated total density vector tr_rho_(r)
      arma::vec tr_rho_;

      /// The on-top pair-density vector pi_(r)
      arma::vec pi_;

      /// The R(r) factor
      arma::vec R_;
};

}
#endif // MCPDFT_H
