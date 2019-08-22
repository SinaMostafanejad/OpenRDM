#include "mcpdft.h"
#include "energy.h"
#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>

using namespace mcpdft;

int main() {
   
    arma::mat D1a, D1b, D2ab; 
    D1a.zeros(); D1b.zeros(); D2ab.zeros();

    double e = mcpdft_energy(D1a,D1b,D2ab);


    printf("==========================================\n");
    printf("   MCPDFT energy =  %-20.15lf\n",   e);
    printf("==========================================\n");
    return 0;
}
