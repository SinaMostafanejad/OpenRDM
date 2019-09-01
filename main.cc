#include "mcpdft.h"
#include "energy.h"
#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>

using namespace mcpdft;

int main() {
  
    MCPDFT* mc = new MCPDFT;

    // fetching the number of basis functions
    int nbfs;
    nbfs = mc->get_nbfs();

    /* building alpha and beta 1-electron reduced
     * density matrices (1RDMs)
     */
    // mc->build_opdm();
    arma::mat D1a(mc->get_D1a());
    arma::mat D1b(mc->get_D1b());

    /* building the alpha-beta block of the 2-electron reduced
     * density matrix (2RDM)
     */
    mc->build_tpdm();
    arma::mat D2ab(mc->get_D2ab());

    mc->build_rho();

    double e = mcpdft_energy(mc,D1a,D1b,D2ab);

    // getting the value of the reference energy
    double eref = mc->get_eref();
    // printf("eref = %-20.15lf\n",eref);

    printf("=================================================\n");
    printf("   Reference energy      =  %-20.15lf\n",  eref);
    printf("   MCPDFT energy         =  %-20.15lf\n",     e);
    printf("   E(MCPDFT) - E(Ref)    =  %-20.2le\n", e-eref);
    printf("=================================================\n");

    delete mc;

    return 0;
}
