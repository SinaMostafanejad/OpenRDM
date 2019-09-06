#include "energy.h"
#include "mcpdft.h"
#include "functional.h"
#include <armadillo>


namespace mcpdft {

   // calculates the MCPDFT energy
   double mcpdft_energy(const MCPDFT* mc,
                        const arma::mat &D1a,
                        const arma::mat &D1b,
                        const arma::mat &D2ab) {

      double tot_energy = 0.0;
      double ref_energy = 0.0;

      // fetching the number of grid points
      size_t npts;
      npts = mc->get_npts();
      // printf("npts = %d\n",(int)npts);

      // fetching the number of basis functions
      int nbfs;
      nbfs = mc->get_nbfs();
      // printf("npts = %d\n",(int)npts);

      // fetching the weights and the coordinates of grids
      arma::vec W(mc->get_w());
      arma::vec X(mc->get_x());
      arma::vec Y(mc->get_y());
      arma::vec Z(mc->get_z());
      // W.print("W = ");
      // X.print("X = ");
      // Y.print("Y = ");
      // Z.print("Z = ");
     
      // getting the values of the basis functions on grids
      arma::mat Phi(mc->get_phi());
      // Phi.print("Phi = ");

      // getting the value of the reference energy
      double eref = mc->get_eref();
      // printf("eref = %-20.15lf\n",eref);

      // getting the value of the classical energy
      double eclass = mc->get_eclass();
      // printf("eclass = %-20.15lf\n",eclass);

      arma::vec rhoa(mc->get_rhoa());
      arma::vec rhob(mc->get_rhob());

      Functional* func = new Functional;

      double Ex = 0.0;
      double Ec = 0.0;
      Ex = func->EX_LSDA(mc, rhoa, rhob);
      Ec = func->EC_VWN3(mc, rhoa, rhob);

      printf("------------------------------------------\n");
      printf("   Classical energy = %-20.15lf\n", eclass);
      printf("   Ex               = %-20.15lf\n", Ex);
      printf("   Ec               = %-20.15lf\n", Ec);
      printf("------------------------------------------\n\n");

      // getting the AO->MO transformation matrix C
      // arma::mat cmat(mc->get_cmat());
      // cmat.print("Cmat = ");

      ref_energy += eref;
      tot_energy += eclass;
      tot_energy += Ex;
      tot_energy += Ec;

      delete func;

      return tot_energy;
   }

}
