#include <armadillo>
#include <sys/sysinfo.h>
#include "energy.h"
#include "mcpdft.h"
#include "libMem.h"
#include "openrdmConfig.h"
#include <string>

#ifdef WITH_LIBXC
   #include <xc.h>
#else
   #include "functional.h"
#endif

namespace mcpdft {

   // calculates the MCPDFT energy
   double mcpdft_energy(MCPDFT *mc,
		        std::string &functional,
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
      delete libmem;

      double tot_energy = 0.0;

      // getting the value of the classical energy
      // double eclass = mc->get_eclass();
      // double eclass = -2.47443074 + 1.34513977 + 0.70556961;

      /* calculate core (kinetic + nuclear attraction) energy */
      arma::mat Hcore(mc->get_hcore());
      arma::mat D1(D1a+D1b);
      double e1el = mc->core_energy(D1,Hcore);
      //printf("e1el = %-20.15lf\n",e1el);

      /* get nuclear repulsion energy */
      double enuc = mc->get_enuc();
      // printf("enuc = %-20.15lf\n",enuc);

      /* classical Hartree electronic repulsion energy */
      arma::mat Ja(mc->get_ja());
      arma::mat Jb(mc->get_jb());
      double eHartree = mc->Hartree_energy(D1a,D1b,Ja,Jb);
      // printf("eHartree = %-20.15lf\n",eHartree);

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
      arma::vec tr_sigma_aa(npts, arma::fill::zeros);
      arma::vec tr_sigma_ab(npts, arma::fill::zeros);
      arma::vec tr_sigma_bb(npts, arma::fill::zeros);
      if(mc->is_gga()) {
         tr_sigma_aa = mc->get_tr_sigma_aa();
         tr_sigma_ab = mc->get_tr_sigma_ab();
         tr_sigma_bb = mc->get_tr_sigma_bb();
      }
#ifdef WITH_LIBXC
      arma::vec tr_rho(tr_rhoa);
      arma::vec ex(npts);
      arma::vec ec(npts);
      arma::vec W(mc->get_w());
      tr_rho = tr_rho + tr_rhob;
      double * rhop = tr_rho.memptr();
      double * exp = ex.memptr();
      double * ecp = ec.memptr();

      double Ex = 0.0;
      double Ec = 0.0;

      xc_func_type func_x;
      if(xc_func_init(&func_x,XC_LDA_X, XC_UNPOLARIZED) != 0){
        xc_lda_exc(&func_x, npts, rhop, exp);
      }

      double ans = 0.0;
      for (int p = 0; p < npts; p++) {
          if (tr_rho(p) > 1.0e-20)
             ans += W(p) * ex(p) * tr_rho(p);
      }
      Ex=ans;

      xc_func_type func_c;
      if(xc_func_init(&func_c,XC_LDA_C_VWN_RPA, XC_UNPOLARIZED) != 0){
        xc_lda_exc(&func_c, npts, rhop, ecp);
      }

      ans = 0.0;
      for (int p = 0; p < npts; p++) {
          if (tr_rho(p) > 1.0e-20)
             ans += W(p) * ec(p) * tr_rho(p);
      }
      Ec=ans;
      xc_func_end(&func_x);
      xc_func_end(&func_c);
#else
      Functional* func = new Functional;

      double Ex = 0.0;
      double Ec = 0.0;
      if (functional == "SVWN") {
          Ex = func->EX_LSDA(mc, tr_rhoa, tr_rhob);
          Ec = func->EC_VWN3(mc, tr_rhoa, tr_rhob);
      }else{
          Ex = func->EX_PBE(mc, tr_rhoa, tr_rhob, tr_sigma_aa, tr_sigma_bb);
          Ec = func->EC_PBE(mc, tr_rhoa, tr_rhob, tr_sigma_aa, tr_sigma_ab, tr_sigma_bb);
      }

      delete func;
#endif

      tot_energy += enuc;
      tot_energy += e1el;
      tot_energy += eHartree;
      tot_energy += Ex;
      tot_energy += Ec;

      printf("-----------------------------------------------------\n");
      printf("   Nuclear repulsion energy = %20.12lf\n", enuc);
      printf("   One-electron energy      = %20.12lf\n", e1el);
      printf("   Classical Hartree energy = %20.12lf\n", eHartree);
      printf("   Ex                       = %20.12lf\n", Ex);
      printf("   Ec                       = %20.12lf\n", Ec);
      printf("-----------------------------------------------------\n");
      printf("   MCPDFT energy            = %20.12lf\n", tot_energy);
      printf("-----------------------------------------------------\n\n");

      return tot_energy;
   }

}
