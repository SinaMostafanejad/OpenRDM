#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <sys/sysinfo.h>
#include "mcpdft.h"
#include "energy.h"
#include "libMem.h"

using namespace mcpdft;

int main(int argc, char *argv[]) {
    if ( argc < 2 ) { 

       std::cout << "An argument is needed!" << std::endl;
       printf("Usage: %s <test_case>\n", argv[0]);
       return 1;
    }

    // Query memory information from linux
    struct sysinfo info;
    sysinfo(&info);
 
    // Calculating the amount of available memory
    LibMem *libmem;
    libmem = new LibMem();
    libmem->query_system_memory(&info);

    // MCPDFT* mc = new MCPDFT(test_case);
    MCPDFT *mc;
    mc = new MCPDFT(argv[1]);

    // getting the value of the reference energy
    double eref = mc->get_eref();
    // printf("eref = %-20.15lf\n",eref);

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

    // calculating the MCPDFT energy correction
    double e = mcpdft_energy(mc,D1a,D1b,D2ab);

    printf("=================================================\n");
    printf("   Reference energy      =  %-20.12lf\n",  eref);
    printf("   MCPDFT energy         =  %-20.12lf\n",     e);
    printf("   E(MCPDFT) - E(Ref)    =  %-20.2le\n", e-eref);
    printf("=================================================\n");

    delete libmem;
    delete mc;

    return 0;
}
