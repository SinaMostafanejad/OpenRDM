#include <armadillo>
#include <sys/sysinfo.h>
//#include <xc.h>
#include "energy.h"
#include "mcpdft.h"
#include "functional.h"
#include "libMem.h"

namespace mcpdft {

   // calculates the MCPDFT energy
   double mcpdft_energy(MCPDFT *mc,
                        const arma::mat &D1a,
                        const arma::mat &D1b,
                        const arma::mat &D2ab) {

      // Query memory information from linux
      struct sysinfo info;
      sysinfo(&info);

      // Calculating the amount of available memory
      LibMem *libmem;
      libmem = new LibMem();
      libmem->query_system_memory(&info);

      double tot_energy = 0.0;

      // getting the value of the classical energy
      double eclass = mc->get_eclass();
      // printf("eclass = %-20.15lf\n",eclass);

      // building the one electron densities rho_a(r) and rho_b(r)
      mc->build_rho();

      // building the on-top pair density pi(r,r)
      mc->build_pi(D2ab);

      // building the R(r) factor for density translation
      mc->build_R();

      // translate the one-electron densities
      mc->translate();

      size_t npts = mc->get_npts();
      arma::vec tr_rhoa(mc->get_tr_rhoa());
      arma::vec tr_rhob(mc->get_tr_rhob());
      //arma::vec tr_rho(tr_rhoa);
      //arma::vec ex(npts);
      //arma::vec ec(npts);
      //arma::vec W(mc->get_w());
      //tr_rho = tr_rho + tr_rhob;
      //double * rhop = tr_rho.memptr();
      //double * exp = ex.memptr();
      //double * ecp = ec.memptr();

      Functional* func = new Functional;

      double Ex = 0.0;
      double Ec = 0.0;
      Ex = func->EX_LSDA(mc, tr_rhoa, tr_rhob);
      Ec = func->EC_VWN3(mc, tr_rhoa, tr_rhob);
      // xc_func_type func_x;
      // if(xc_func_init(&func_x,XC_LDA_X, XC_UNPOLARIZED) != 0){
      //   xc_lda_exc(&func_x, npts, rhop, exp);
      // }

      // double ans = 0.0;
      // for (int p = 0; p < npts; p++) {
      //     if (tr_rho(p) > 1.0e-20)
      //        ans += W(p) * ex(p) * tr_rho(p);
      // }
      // Ex=ans;

      // xc_func_type func_c;
      // if(xc_func_init(&func_c,XC_LDA_C_VWN_RPA, XC_UNPOLARIZED) != 0){
      //   xc_lda_exc(&func_c, npts, rhop, ecp);
      // }

      // ans = 0.0;
      // for (int p = 0; p < npts; p++) {
      //     if (tr_rho(p) > 1.0e-20)
      //        ans += W(p) * ec(p) * tr_rho(p);
      // }
      // Ec=ans;
      // xc_func_end(&func_x);
      // xc_func_end(&func_c);

      printf("------------------------------------------\n");
      printf("   Classical energy = %-20.12lf\n", eclass);
      printf("   Ex               = %-20.12lf\n", Ex);
      printf("   Ec               = %-20.12lf\n", Ec);
      printf("------------------------------------------\n\n");

      tot_energy += eclass;
      tot_energy += Ex;
      tot_energy += Ec;

      delete func;
      delete libmem;

      return tot_energy;
   }

}
