#include <armadillo>
#include <sys/sysinfo.h>
#include "energy.h"
#include "mcpdft.h"
#include "libMem.h"
#include "openrdmConfig.h"
#include <string>
#include "HDF5_Read_Contiguous.h"
#include "HDF5_Write_Contiguous.h"
#include "HDF5Factory.h"
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

      size_t nbfs = mc->get_nbfs();
      size_t nbfs2 = nbfs * nbfs;
      arma::mat d1a(nbfs, nbfs, arma::fill::zeros);
      arma::mat d1b(nbfs, nbfs, arma::fill::zeros);
      arma::mat d2ab(nbfs2, nbfs2, arma::fill::zeros);
      IOFactory* iof;
      IRead* ird;
      IWrite* iwt;

      iof = new HDF5Factory;
      iwt = iof->create_IWrite();
      iwt->write_rdms(D1a,D1b,D2ab);
      ird = iof->create_IRead();
      ird->read_rdms(d1a,d1b,d2ab);
      delete iof;
      // IWrite* h5w = iof->create_IWrite();
      // h5w->write_opdm(D1a,D1b);
      // DiskRW dskrw;
      // dskrw.write_opdm(D1a,D1b);
      // dskrw.write_tpdm(D2ab);
      // size_t nbfs = mc->get_nbfs();
      // size_t nbfs2 = nbfs * nbfs;
      // try{
      //    arma::mat d1a(nbfs, nbfs, arma::fill::zeros);
      //    arma::mat d1b(nbfs, nbfs, arma::fill::zeros);
      //    arma::mat d2ab(nbfs2, nbfs2, arma::fill::zeros);
      //    dskrw.read_opdm(d1a,d1b);
      //    // d1a.print("D1a =");
      //    // d1b.print("D1b =");
      //    dskrw.read_tpdm(d2ab);
      //    // d2ab.print("D2ab =");
      // } catch(const char* err_msg) {
      //    printf("%s\n",err_msg);
      // }


      printf("------------------------------------------\n");
      printf("   Classical energy = %-20.12lf\n", eclass);
      printf("   Ex               = %-20.12lf\n", Ex);
      printf("   Ec               = %-20.12lf\n", Ec);
      printf("------------------------------------------\n\n");

      tot_energy += eclass;
      tot_energy += Ex;
      tot_energy += Ec;

      delete libmem;

      return tot_energy;
   }

}
