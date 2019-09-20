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

      // getting the value of the reference energy
      double eref = mc->get_eref();
      // printf("eref = %-20.15lf\n",eref);

      // getting the value of the classical energy
      double eclass = mc->get_eclass();
      // printf("eclass = %-20.15lf\n",eclass);

      arma::vec tr_rhoa(mc->get_tr_rhoa());
      arma::vec tr_rhob(mc->get_tr_rhob());

      Functional* func = new Functional;

      double Ex = 0.0;
      double Ec = 0.0;
      Ex = func->EX_LSDA(mc, tr_rhoa, tr_rhob);
      Ec = func->EC_VWN3(mc, tr_rhoa, tr_rhob);

      printf("------------------------------------------\n");
      printf("   Classical energy = %-20.12lf\n", eclass);
      printf("   Ex               = %-20.12lf\n", Ex);
      printf("   Ec               = %-20.12lf\n", Ec);
      printf("------------------------------------------\n\n");

      ref_energy += eref;
      tot_energy += eclass;
      tot_energy += Ex;
      tot_energy += Ec;

      delete func;

      return tot_energy;
   }

}
