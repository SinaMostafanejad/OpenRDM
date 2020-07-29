#include <stdio.h>
#include <armadillo>
#include "mcpdft.h"
#include "energy.h"
#include "erefs.h"

using namespace mcpdft;

int main(int argc, char *argv[]) {
    if ( argc < 2 ) { 
       std::cout << "An argument is needed!" << std::endl;
       printf("Usage: %s <test_case> <functional>\n", argv[0]);
       return 1;
    }

    // MCPDFT* mc = new MCPDFT(test_case);
    MCPDFT *mc;
    // mc = new MCPDFT(argv[1]);
    mc = new MCPDFT();

    // getting the value of the reference energy
    std::string functional(argv[2]);
    double eref{0.0};
    if (functional == "SVWN") {
       eref = 0.0; 
    }else{
       eref = -1.156522359214;
    }

    /* fetching alpha and beta 1-electron reduced
     * density matrices (1-RDMs)
     */
    arma::mat D1a(mc->get_D1a());
    arma::mat D1b(mc->get_D1b());

    /* fetching the alpha-beta block of the 2-electron reduced
     * density matrix (2-RDM)
     */
    arma::mat D2ab(mc->get_D2ab());

    // calculating the MCPDFT energy correction
    double e =  mcpdft_energy(mc,functional,D1a,D1b,D2ab);

    printf("=================================================\n");
    printf("   Reference energy      =  %-20.12lf\n",  eref);
    printf("   MCPDFT energy         =  %-20.12lf\n",     e);
    printf("   E(MCPDFT) - E(Ref)    =  %-20.2le\n", e-eref);
    printf("=================================================\n\n");

    delete mc;

    return 0;
}
