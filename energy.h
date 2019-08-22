#ifndef ENERGY_H
#define ENERGY_H

#include <armadillo>
#include <iostream>
#include <fstream>

#include "mcpdft.h"

namespace mcpdft {

   /// calculates the MCPDFT energy
   double mcpdft_energy(const MCPDFT* mc,
                        const arma::mat &D1a,
                        const arma::mat &D1b,
                        const arma::mat &D2ab);
}

#endif // ENERGY_H
