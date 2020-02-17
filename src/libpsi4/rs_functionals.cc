/*
 * @BEGIN LICENSE
 *
 * mydft by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>
#include "psi4/libpsi4util/libpsi4util.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// for grid
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

// jk object
#include "psi4/libfock/jk.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

#include "mcpdft_solver.h"

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ The four major conventions regarding to the principal variables used for                                              @
// @ building exchange-correlation (XC) functionals are as follows:                                                        @
// @                                                                                                                       @
// @ Convention (I)   spin densities rho_a, rho_b and their gradients:                   EXC[rho_a, rho_b, rho_a', rho_b'] @
// @ Convention (II)  total density, spin-magnetization density m and their gradients:   EXC[rho, m, rho', m']             @
// @ Convention (III) total density, spin polarization and their gradients:              EXC[rho, zeta, rho', zeta']       @
// @ Convention (IV)  singlet and triplet charge densities                                                                 @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

namespace psi{ namespace mcpdft {

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Exchange functionals ++++++++++++++++++
    // wPBE                                                     +
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MCPDFTSolver::EX_wPBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                             std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB) {

    const double OMEGA = options_.get_double("MCPDFT_OMEGA");

    const double A_bar = 0.757211;
    const double B = -0.106364;
    const double C = -0.118649;
    const double D = 0.609650;
    const double E = -0.0477963;

    const double a2 = 0.0159941;
    const double a3 = 0.0852995;
    const double a4 = -0.160368;
    const double a5 = 0.152645;
    const double a6 = -0.0971263;
    const double a7 = 0.0422061;
    const double b1 = 5.33319;
    const double b2 = -12.4780;
    const double b3 = 11.0988;
    const double b4 = -5.11013;
    const double b5 = 1.71468;
    const double b6 = -0.610380;
    const double b7 = 0.307555;
    const double b8 = -0.0770547;
    const double b9 = 0.0334840;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    auto kF = [](double RHO) -> double {

              double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
              return dum;
    };

    auto S = [=](double RHO, double SIGMA) -> double {

             double temp = sqrt(SIGMA) / (2.0 * kF(RHO) * RHO);
             return temp;
    };

    auto H = [=](double RHO, double SIGMA) -> double {

             double s1 = S(RHO,SIGMA);
             double s2 = s1 * s1;
             double s3 = s2 * s1;
             double s4 = s3 * s1;
             double s5 = s4 * s1;
             double s6 = s5 * s1;
             double s7 = s6 * s1;
             double s8 = s7 * s1;
             double s9 = s8 * s1;

             double numerator   = a2 * s2 + a3 * s3 + a4 * s4 + a5 * s5 + a6 * s6 + a7 * s7;
             double denominator = 1.0 + b1 * s1 + b2 * s2 + b3 * s3 + b4 * s4 + b5 * s5 + b6 * s6 + b7 * s7 + b8 * s8 + b9 * s9;
             double dum = numerator/denominator;

             return dum;
    };

    auto Nu = [=](double RHO) -> double {

              double dum = OMEGA / kF(RHO);
              return dum;
    };

    auto ZETA = [=](double RHO, double SIGMA) -> double {

             double dum = pow(S(RHO,SIGMA), 2.0) * H(RHO, SIGMA);
             return dum;

    };

    auto ETA = [=](double RHO, double SIGMA) -> double {

             double temp = A_bar + ZETA(RHO, SIGMA);
             return temp;

    };

    auto LAMBDA = [=](double RHO, double SIGMA) -> double {

             double temp = D + ZETA(RHO, SIGMA);
             return temp;

    };

    auto Chi = [=](double RHO, double SIGMA) -> double {

              double dum = Nu(RHO) / sqrt(LAMBDA(RHO,SIGMA) + pow(Nu(RHO),2.0));
              return dum;
    };

    auto BG_LAMBDA = [=](double RHO, double SIGMA) -> double {

             double temp = LAMBDA(RHO,SIGMA) / (1.0 - Chi(RHO,SIGMA));
             return temp;

    };

    auto F_bar = [=](double RHO, double SIGMA) -> double {

                 double s0 = 2.0;
                 double s1 = S(RHO,SIGMA);
                 double s2 = s1 * s1;
                 double ze_v = ZETA(RHO,SIGMA);

                 double dum = 1.0 - (1.0/(27.0*C)) * (s2/(1.0 + (s2/pow(s0,2.0)))) - (1.0/(2.0*C)) * ze_v;

                 return dum;

    };

    auto G_bar = [=](double RHO, double SIGMA) -> double {

                 double ze_v   = ZETA(RHO,SIGMA);
                 double et_v   = ETA(RHO,SIGMA);
                 double lam_v  = LAMBDA(RHO,SIGMA);
                 double lam_v2 = lam_v * lam_v;
                 double lam_v3 = lam_v2 * lam_v;
                 double lam_v72 = pow(lam_v,7.0/2.0);
                 double sq_ze = sqrt(ze_v);
                 double sq_et = sqrt(et_v);

                 double dum = -(2.0/5.0) * C * F_bar(RHO,SIGMA) * lam_v - (4.0/15.0) * B * lam_v2 - (6.0/5.0) * A_bar * lam_v3
                            - (4.0/5.0) * sqrt(M_PI) * lam_v72 - (12.0/5.0) * lam_v72 * (sq_ze - sq_et);

                 dum *= (1.0/E);

                 return dum;

    };

    auto eX = [=](double RHO) -> double {

              double temp = -(3.0 * kF(RHO)) / (4.0 * M_PI);
              return temp;
    };

    auto FX = [=](double RHO, double SIGMA) -> double {

              double chi_v  = Chi(RHO,SIGMA);
              double chi_v2 = chi_v * chi_v;
              double chi_v3 = chi_v2 * chi_v;
              double chi_v5 = chi_v3 * chi_v2;
              double nu_v   = Nu(RHO);
              double ze_v   = ZETA(RHO,SIGMA);
              double et_v   = ETA(RHO,SIGMA);
              double bg_lam_v  = BG_LAMBDA(RHO,SIGMA);
              double bg_lam_v2 = bg_lam_v * bg_lam_v;
              double bg_lam_v3 = bg_lam_v2 * bg_lam_v;
              double lam_v  = LAMBDA(RHO,SIGMA);
              double lam_v2 = lam_v * lam_v;
              double lam_v3 = lam_v2 * lam_v;
              double sq_ze_nu = sqrt(ze_v  + nu_v*nu_v);
              double sq_et_nu = sqrt(et_v  + nu_v*nu_v);
              double sq_la_nu = sqrt(lam_v + nu_v*nu_v);

              double temp = A_bar * ( ze_v/( (sq_ze_nu + sq_et_nu) * (sq_ze_nu + nu_v) ) + et_v/( (sq_ze_nu + sq_et_nu) * (sq_et_nu + nu_v) ) )
                          - (4.0/9.0) * (B/bg_lam_v) - (4.0/9.0) * (C * F_bar(RHO,SIGMA) / bg_lam_v2) * (1.0 + (1.0/2.0) * chi_v )
                          - (8.0/9.0) * (E * G_bar(RHO,SIGMA) / bg_lam_v3) * (1.0 + (9.0/8.0) * chi_v + (3.0/8.0) * chi_v2)
                          + 2.0 * ze_v * log( 1.0 - (lam_v - ze_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_ze_nu) ) )
                          - 2.0 * et_v * log( 1.0 - (lam_v - et_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_et_nu) ) );

              return temp;
    };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;

        double sigmaaa = sigma_aap[p];
        double sigmabb = sigma_bbp[p];
        double sigma = 0.0;

        double tol = 1.0e-20;
        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;

              double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
              exc += rho * zk * grid_w_->pointer()[p];

           }else if ( rhob < tol ){

                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;

                    double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
                    exc += rho * zk * grid_w_->pointer()[p];
           }else {

                 double zka = rhoa * eX(2.0 * rhoa) * FX(2.0 * rhoa, 4.0 * sigmaaa);
                 double zkb = rhob * eX(2.0 * rhob) * FX(2.0 * rhob, 4.0 * sigmabb);
                 double zk = zka + zkb;
                 exc += zk * grid_w_->pointer()[p];
           }
        }else{
                //double zk = 0.0;        
                exc += 0.0;
             }
    }
    return exc;
}

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++ Correlation Functionals +++++++++++++++++
    // 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


}} // End namespaces
