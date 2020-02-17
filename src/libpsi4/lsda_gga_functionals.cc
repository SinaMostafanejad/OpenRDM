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

namespace psi{ namespace mcpdft {

    //###########################################################
    //# The exchange and correlation functional implementations #
    //###########################################################

double MCPDFTSolver::Gfunction(double r, double A, double a1, double b1, double b2, double b3, double b4, double p) {

    double G = -2.0 * A * (1.0 + a1 * r) * log( 1.0 + pow( 2.0 * A * ( b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0 ) )  ,-1.0 ) );

    return G;

}
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Exchange functionals ++++++++++++++++++
    // LDA/LSDA, B86_MGC, B88, PBE                              +
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


double MCPDFTSolver::EX_LDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);
    
    double * rho_ap = rho_a->pointer();
    double * rho_bp = rho_b->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_ap[p] + rho_bp[p];

        exc += -Cx * pow( rho, 4.0/3.0) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_LSDA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double exa = pow(2.0,1.0/3.0) * Cx * pow( rho_ap[p], 4.0/3.0) ;
        double exb = pow(2.0,1.0/3.0) * Cx * pow( rho_bp[p], 4.0/3.0) ;
        double ex_LSDA = exa + exb;
        exc += - ex_LSDA * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_LSDA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * zeta_p = ZETA->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
   
        // build f(zeta) weight factor where f(0) = 0 and f(1) = 1
        double fZet =  ( pow( (1.0 + zeta_p[p]) ,4.0/3.0 ) + pow( (1.0 - zeta_p[p]) ,4.0/3.0) - 2.0) / ( 2.0 * pow(2,1.0/3.0) - 2.0 );

        double ex0 = Cx * pow( (rho_ap[p] + rho_bp[p]), 1.0/3.0) ;
        double ex1 = pow(2.0,1.0/3.0) * ex0;
        double ex_LSDA = ex0 + (ex1 - ex0) * fZet;
        exc += - ex_LSDA * (rho_ap[p] + rho_bp[p]) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_B86_MGC(){
    
    const double Cx = 0.73855876638202240586; 
    const double beta = 0.00375;
    const double c = pow(2.0,1.0/3.0) * Cx;
   
    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();
 
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        double Xa_2 = Xa * Xa;
        double Xb_2 = Xb * Xb;
         
        exc += ( -c * rhoa_43 - (beta * Xa_2 * rhoa_43) / pow(1.0 + 0.007 * Xa_2,4.0/5.0) ) * grid_w_->pointer()[p]; 
        exc += ( -c * rhob_43 - (beta * Xb_2 * rhob_43) / pow(1.0 + 0.007 * Xb_2,4.0/5.0) ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_B88(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                            std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){
    
    const double Cx = 0.73855876638202240586; 
    const double beta = 0.0042;
    const double c = pow(2.0,1.0/3.0) * Cx;
 
    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();
 
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        double Xa_2 = Xa * Xa;
        double Xb_2 = Xb * Xb;

        exc += -rhoa_43 * ( c + (beta * Xa_2) / (1.0 + 6.0 * beta * Xa * asinh(Xa)) ) * grid_w_->pointer()[p]; 
        exc += -rhob_43 * ( c + (beta * Xb_2) / (1.0 + 6.0 * beta * Xb * asinh(Xb)) ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_B88_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                           std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){

    double tol = 1.0e-20;

    const double Cx = 0.73855876638202240586;
    const double beta = 0.0042;
    const double c = pow(2.0,1.0/3.0) * Cx;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double sigmaaa = sigma_aap[p];
        double sigmabb = sigma_bbp[p];
        double rhoa_43 = 0.0;
        double rhob_43 = 0.0;
        double Xa   = 0.0;
        double Xb   = 0.0;
        double Xa_2 = 0.0;
        double Xb_2 = 0.0;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rhoa = 0.0;
              sigmabb = std::max(0.0,sigma_bbp[p]);
              rhob_43 = pow( rhob, 4.0/3.0);
              Xb = sqrt(sigma_bbp[p]) / rhob_43;
              Xb_2 = Xb * Xb;

           }else if ( rhob < tol ){

                    rhob = 0.0;
                    sigmaaa = std::max(0.0,sigma_aap[p]);
                    rhoa_43 = pow( rhoa, 4.0/3.0);
                    Xa = sqrt(sigma_aap[p]) / rhoa_43;
                    Xa_2 = Xa * Xa;
           }else{

                sigmaaa = std::max(0.0,sigma_aap[p]);
                sigmabb = std::max(0.0,sigma_bbp[p]);
                rhoa_43 = pow( rhoa, 4.0/3.0);
                rhob_43 = pow( rhob, 4.0/3.0);
                Xa = sqrt(sigma_aap[p]) / rhoa_43;
                Xb = sqrt(sigma_bbp[p]) / rhob_43;
                Xa_2 = Xa * Xa;
                Xb_2 = Xb * Xb;

           }
           exc += -rhoa_43 * ( c + (beta * Xa_2) / (1.0 + 6.0 * beta * Xa * asinh(Xa)) ) * grid_w_->pointer()[p];
           exc += -rhob_43 * ( c + (beta * Xb_2) / (1.0 + 6.0 * beta * Xb * asinh(Xb)) ) * grid_w_->pointer()[p];

        }else{

             exc += 0.0;
        }
    }
    return exc;
}

double MCPDFTSolver::EX_PBE(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                            std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB) {

    const double delta = 0.06672455060314922;
    const double MU = (1.0/3.0) * delta * M_PI * M_PI;
    const double KAPPA = 0.804;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
   
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        
        double Sa = (Xa * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
        double Sb = (Xb * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
        
        double Fsa = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sa,2.0)) / KAPPA ), -1.0 );
        double Fsb = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sb,2.0)) / KAPPA ), -1.0 );
       
        auto E = [](double rhos, double Fss) -> double{

                 double temp = -0.75 * pow(3.0, 1.0/3.0) * pow(M_PI, 2.0/3.0) * pow(rhos,4.0/3.0) * Fss / M_PI;
                 return temp;
        };
 
        double EX_GGAa = 0.5 * E(2.0*rhoa,Fsa);
        double EX_GGAb = 0.5 * E(2.0*rhob,Fsb);
       
        exc += ( EX_GGAa + EX_GGAb ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_PBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                           std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB) {


    const double delta = 0.06672455060314922;
    const double MU = (1.0/3.0) * delta * M_PI * M_PI;
    const double KAPPA = ( (options_.get_str("MCPDFT_FUNCTIONAL") == "REVPBE") ? 1.245 : 0.804 );

    double tol = 1.0e-20;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    auto kF = [](double RHO) -> double {

              double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
              return dum;
    };

    auto eX = [=](double RHO) -> double {

              double temp = -(3.0 * kF(RHO)) / (4.0 * M_PI);
              return temp;
    };

    auto FX = [=](double SIGMA) -> double {

              double temp = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(SIGMA,2.0)) / KAPPA ), -1.0 );
              return temp;
    };


    auto S = [=](double RHO, double SIGMA) -> double {

             double temp = sqrt(SIGMA) / (2.0 * kF(RHO) * RHO);
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

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;

              double zk = eX(2.0 * rho) * FX(S(2.0 * rho, 4.0 * sigma));
              exc += rho * zk * grid_w_->pointer()[p];

           }else if ( rhob < tol ){

                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;

                    double zk = eX(2.0 * rho) * FX(S(2.0 * rho, 4.0 * sigma));
                    exc += rho * zk * grid_w_->pointer()[p];
           }else {

                 double zka = rhoa * eX(2.0 * rhoa) * FX(S(2.0 * rhoa, 4.0 * sigmaaa));
                 double zkb = rhob * eX(2.0 * rhob) * FX(S(2.0 * rhob, 4.0 * sigmabb));
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

//
// 
//     const double alpha = (2.0/3.0);      // Slater value
//     const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);
//     const double MU = 0.2195149727645171;
//     const double KAPPA = 0.804;
// 
//     double * rho_ap = rho_a_->pointer();
//     double * rho_bp = rho_b_->pointer();
// 
//     double * rho_a_xp = rho_a_x_->pointer();
//     double * rho_b_xp = rho_b_x_->pointer();
// 
//     double * rho_a_yp = rho_a_y_->pointer();
//     double * rho_b_yp = rho_b_y_->pointer();
// 
//     double * rho_a_zp = rho_a_z_->pointer();
//     double * rho_b_zp = rho_b_z_->pointer();
// 
//     double * zeta_p = zeta_->pointer();
//     double * sigma_aap = sigma_aa_->pointer();
//     double * sigma_abp = sigma_ab_->pointer();
//     double * sigma_bbp = sigma_bb_->pointer();
// 
//     double exc = 0.0;
//     for (int p = 0; p < phi_points_; p++) {
// 
//         double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
//         // double sig = sigma_aap[p]; 
//         double rhoa = rho_ap[p];
//         double rhob = rho_bp[p];
//         double rho = rhoa + rhob;
//         // local fermi wave vector
//         double kf = pow( ( 3.0 * pow(M_PI,2.0) * rho), 1.0/3.0);
//         // double kfa = pow( ( 3.0 * pow(M_PI,2.0) * rhoa ) , 1.0/3.0);
//         // double kfb = pow( ( 3.0 * pow(M_PI,2.0) * rhob ) , 1.0/3.0);
// 
//         // double EXa = -(3.0 * kfa) / (4.0 * M_PI); 
//         // double EXb = -(3.0 * kfb) / (4.0 * M_PI);
//         // double absDelRho = sqrt( ( pow( (rho_a_xp[p] + rho_b_xp[p]) ,2.0) + pow( (rho_a_yp[p] + rho_b_yp[p]) ,2.0) + pow( (rho_a_zp[p] + rho_b_zp[p]) ,2.0) ) );
//         double absDelRho = sqrt( sig );
// 
//         // double Sa = sqrt(4.0 * grada_2) / (2.0 * kfa * rhoa); 
//         // double Sb = sqrt(4.0 * gradb_2) / (2.0 * kfb * rhob);
// 
//         double s = absDelRho / ( 2.0 * kf * rho);
//         // double sa = absDelRho / ( 2.0 * kfa *2.0* rho_ap[p] );
//         // double sb = absDelRho / ( 2.0 * kfb *2.0* rho_bp[p] );
// 
//         double Fs = pow( (1.0 + 1.296 * pow(s,2.0) + 14.0 * pow(s,4.0) + 0.2 * pow(s,6.0) ) , 1.0/15.0 ) ;
//         // double Fsa = pow( (1.0 + 1.296 * pow(sa,2.0) + 14.0 * pow(sa,4.0) + 0.2 * pow(sa,6.0) ) , 1.0/15.0 );
//         // double Fsb = pow( (1.0 + 1.296 * pow(sb,2.0) + 14.0 * pow(sb,4.0) + 0.2 * pow(sb,6.0) ) , 1.0/15.0 );
//         // double Fs = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(s,2.0)) / KAPPA ), -1.0 );
//         // double Fsa = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sa,2.0)) / KAPPA ), -1.0 );
//         // double Fsb = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sb,2.0)) / KAPPA ), -1.0 );
// 
//         // double EX_GGAa = rhoa * EXa * Fsa;
//         // double EX_GGAb = rhob * EXb * Fsb;
//         exc += -Cx * pow( rho, 4.0/3.0) * Fs * grid_w_->pointer()[p];
//         // exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * Fs * grid_w_->pointer()[p]; 
//         // exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) * Fsa + pow( rho_bp[p], 4.0/3.0) * Fsb ) * grid_w_->pointer()[p]; 
//         // exc += -0.5 * ( -3.0 *2.0* kfa/(4*M_PI) * rho_ap[p] * Fsa - 3.0 *2.0* kfb / (4*M_PI) * rho_bp[p] * Fsb ) * grid_w_->pointer()[p]; 
//         // exc += 0.5 * ( EX_GGAa + EX_GGAb ) * grid_w_->pointer()[p]; 
//     }
//     return exc;
// }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++ Correlation Functionals +++++++++++++++++
    // B88_OP, VWN3, LYP, PW92, PBE
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Note: This correlation functional depends on B86MGC exchange functional
// with empirical atomic parameters t and u defined below.
// From T. Tsuneda, J. Chem. Phys. 110, 10664 (1999).

double MCPDFTSolver::EC_B88_OP(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){

   const double beta = 0.0042;
 
   double * rho_ap = RHO_A->pointer();
   double * rho_bp = RHO_B->pointer();
   double * sigma_aap = SIGMA_AA->pointer();
   double * sigma_bbp = SIGMA_BB->pointer();

   double exc = 0.0;
   for (int p = 0; p < phi_points_; p++) {

       double rhoa = rho_ap[p];
       double rhob = rho_bp[p];
       double rhoa_43 = pow( rhoa, 4.0/3.0);
       double rhob_43 = pow( rhob, 4.0/3.0);
       double rhoa_13 = pow( rhoa, 1.0/3.0);
       double rhob_13 = pow( rhob, 1.0/3.0);
       double sigmaaa = sigma_aap[p];
       double sigmabb = sigma_bbp[p];
       double Xa = sqrt(sigmaaa) / rhoa_43;
       double Xb = sqrt(sigmabb) / rhob_43;
       double Xa_2 = Xa * Xa;
       double Xb_2 = Xb * Xb;
      
       auto Ks = [=](double Xs, double Xs_2) -> double{
               
                 double temp3 = 3.0 * pow(3.0/(4.0 * M_PI),1.0/3.0) + 2.0 * (beta * Xs_2) / (1.0 + 6.0 * beta * Xs * asinh(Xs));
                 return temp3; 
       };

       double Ka = Ks(Xa, Xa_2);
       double Kb = Ks(Xb, Xb_2);

       double BETA_ab = 2.3670 * (rhoa_13 * rhob_13 * Ka * Kb) / (rhoa_13 * Ka + rhob_13 * Kb);
       
       exc += -rhoa * rhob * ( (1.5214 * BETA_ab + 0.5764) / ( pow(BETA_ab,4.0) + 1.1284 * pow(BETA_ab,3.0) + 0.3183 * pow(BETA_ab,2.0) ) ) * grid_w_->pointer()[p];

   }
   return exc;
}

// double MCPDFTSolver::EC_B88(){
// 
//    const double Cx = 0.73855876638202240586;
//    const double c = pow(2.0,1.0/3.0) * Cx;
//    const double t = 0.63;
//    const double u = 0.96;
//    const double beta = 0.00375;
//    // const double beta = 0.0042;
//    const double lambda = 0.007;
//  
//    double * rho_ap = rho_a_->pointer();
//    double * rho_bp = rho_b_->pointer();
//    double * sigma_aap = sigma_aa_->pointer();
//    double * sigma_bbp = sigma_bb_->pointer();
//    double * tau_ap = tau_a_->pointer();
//    double * tau_bp = tau_b_->pointer();
// 
//    double exc = 0.0;
//    for (int p = 0; p < phi_points_; p++) {
// 
//        double rhoa = rho_ap[p];
//        double rhob = rho_bp[p];
//        double rho = rhoa + rhob;
//        double rhoa_43 = pow( rhoa, 4.0/3.0);
//        double rhob_43 = pow( rhob, 4.0/3.0);
//        double rhoa_13 = pow( rhoa, 1.0/3.0);
//        double rhob_13 = pow( rhob, 1.0/3.0);
//        double taua = tau_ap[p];
//        double taub = tau_bp[p];
//        double sigmaaa = sigma_aap[p];
//        double sigmabb = sigma_bbp[p];
//        double Xa = sqrt(sigmaaa) / rhoa_43;
//        double Xb = sqrt(sigmabb) / rhob_43;
//        double Xa_2 = Xa * Xa;
//        double Xb_2 = Xb * Xb;
// 
//        auto xyfunc = [=](double rhos_13, double Xs_2) -> double{
//                 
//                      double temp = 0.5 * pow(c * rhos_13 + (beta * Xs_2 * rhos_13) / pow(1.0 + lambda * Xs_2 ,4.0/5.0 ) ,-1.0);
//                      return temp;
//        };
//        
//        auto rfunc = [=](double rhos, double rhos_43, double Xs_2) -> double{
//        
//                 double temp1 =  0.5 * rhos * pow( c * rhos_43 + (beta * Xs_2 * rhos_43)/ pow(1.0 + lambda * Xs_2 ,4.0/5.0) ,-1.0);
//                 return temp1;
//        };
// 
//        auto dfunc = [=](double rhos, double taus, double sigmass) -> double {
// 
//                 double temp2 = taus - 0.25 * sigmass / rhos;
//                 return temp2;
//        };
// 
//        double x = xyfunc(rhoa_13,Xa_2);
//        double y = xyfunc(rhob_13,Xb_2);
//        double q = t * (x + y);
//        double q_2 = q * q;
//        double r_a = rfunc(rhoa, rhoa_43, Xa_2);
//        double r_b = rfunc(rhob, rhob_43, Xb_2);
//        double z_a = 2.0 * u * r_a;
//        double z_b = 2.0 * u * r_b;
//        double d_a = dfunc(rhoa, taua, sigmaaa);
//        double d_b = dfunc(rhob, taub, sigmabb);
//        
//        double f = -0.8 * rhoa * rhob * q_2 * (1.0 - log(1.0 + q) / q);
//        double g_a = -0.01 * rhoa * d_a * pow(z_a,4.0) * ( 1.0 - 2.0 * log(1.0 + 0.5 * z_a) / z_a );
//        double g_b = -0.01 * rhob * d_b * pow(z_b,4.0) * ( 1.0 - 2.0 * log(1.0 + 0.5 * z_b) / z_b );
//       
//        auto Ks = [=](double Xs, double Xs_2) -> double{
//                
//                  double temp3 = 3.0 * pow(3.0/(4.0 * M_PI),1.0/3.0) + 2.0 * (beta * Xs_2) / (1.0 + 6.0 * beta * Xs * asinh(Xs));
//                  return temp3; 
//        };
// 
//        double Ka = Ks(Xa, Xa_2);
//        double Kb = Ks(Xb, Xb_2);
// 
//        exc += ( f + g_a + g_b ) * grid_w_->pointer()[p];
//        double BETA_ab = 2.3670 * (rhoa_13 * rhob_13 * Ka * Kb) / (rhoa_13 * Ka + rhob_13 * Kb);
//        // exc += -rhoa * rhob * ( (1.5214 * BETA_ab + 0.5764) / ( pow(BETA_ab,4.0) + 1.1284 * pow(BETA_ab,3.0) + 0.3183 * pow(BETA_ab,2.0) ) ) * grid_w_->pointer()[p];
// 
//    }
//    return exc;
// }

double MCPDFTSolver::EC_VWN3_RPA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA, std::shared_ptr<Vector> RS){
  
    const double k1 = 0.0310907;
    const double k2 = 0.01554535;
    const double l1 = -0.409286;
    const double l2 = -0.743294;
    const double m1 = 13.0720;
    const double m2 = 20.1231;
    const double n1 = 42.7198;
    const double n2 = 101.578;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    
    double * zeta_p = ZETA->pointer();
    double * rs_p = RS->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = zeta_p[p];
        double rs = rs_p[p];
        double x = sqrt(rs);
        
        double y = (9.0/8.0) * pow(1.0 + zeta, 4.0/3.0) + (9.0/8.0) * pow(1.0 - zeta, 4.0/3.0) - (9.0/4.0);
        double z = (4.0 * y) / (9.0 * pow(2.0,1.0/3.0) - 9.0);

        auto X = [](double i, double c, double d) -> double{
                 
                 double temp = pow(i,2.0) + c * i + d;
                 return temp;
        };
        
        auto Q = [](double c, double d) -> double{
        
                 double temp1 = sqrt( 4 * d - pow(c,2.0) );
                 return temp1;
        };

        auto q = [=](double A, double p, double c, double d) -> double{
            
                 double dum1 = A * ( log( pow(x,2.0) / X(x,c,d) ) + 2.0 * c * atan( Q(c,d)/(2.0*x + c) ) * pow(Q(c,d),-1.0)    
                           - c * p * ( log( pow(x-p,2.0) / X(x,c,d) ) + 2.0 * (c + 2.0 * p) * atan( Q(c,d)/(2.0*x + c) ) * pow(Q(c,d),-1.0) ) * pow(X(p,c,d),-1.0) ); 
                 return dum1;
        };

        double Lambda = q(k1, l1, m1, n1);
        double lambda = q(k2, l2, m2, n2);
       
        double e =  Lambda + z * (lambda - Lambda); 
 
        exc += e * rho * grid_w_->pointer()[p]; 
    }
    return exc;
    
}

// Conv. (III)
double MCPDFTSolver::EC_VWN3_RPA_III(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B){ //, std::shared_ptr<Vector> ZETTA){

    double tol = 1.0e-20;

    const double ecp1 = 0.03109070000;
    const double ecp2 = -0.409286;
    const double ecp3 = 13.0720;
    const double ecp4 = 42.7198;
    const double ecf1 = 0.01554535000;
    const double ecf2 = -0.743294;
    const double ecf3 = 20.1231;
    const double ecf4 = 101.578;
    const double d2Fz = 1.7099209341613656173;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    auto x = [](double RHO) -> double {

             double rs = pow( 3.0 / ( 4.0 * M_PI * RHO ) , 1.0/3.0 );
             double dum = sqrt(rs);
             return dum;
    };

    auto Fz = [](double ZETA) -> double {

              double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
              return dum;
    };

    auto X = [](double i, double c, double d) -> double{

             double temp = pow(i,2.0) + c * i + d;
             return temp;
    };

    auto Q = [](double c, double d) -> double{

             double temp1 = sqrt( 4 * d - pow(c,2.0) );
             return temp1;
    };

    auto q = [=](double RHO, double A, double p, double c, double d) -> double{

             double dum1 = A * ( log( pow(x(RHO),2.0) / X(x(RHO),c,d) ) + 2.0 * c * atan( Q(c,d)/(2.0*x(RHO) + c) ) * pow(Q(c,d),-1.0)
                         - c * p * ( log( pow(x(RHO)-p,2.0) / X(x(RHO),c,d) ) + 2.0 * (c + 2.0 * p) * atan( Q(c,d)/(2.0*x(RHO) + c) )
                         * pow(Q(c,d),-1.0) ) * pow(X(p,c,d),-1.0) );
             return dum1;
    };

    auto EcP = [=](double RHO) -> double {

               double dumm = q(RHO,ecp1,ecp2,ecp3,ecp4);
               return dumm;
    };

    auto EcF = [=](double RHO) -> double {

               double dum = q(RHO,ecf1,ecf2,ecf3,ecf4);
               return dum;
    };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = 0.0;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              zeta = 1.0;

           }else if ( rhob < tol ){

                    rho = rhoa;
                    zeta = 1.0;

           }else {/* if (!(rhoa < tol) && !(rhob < tol) ) */

                 zeta = (rhoa - rhob) / rho;
           }
           double zk = EcP(rho) + Fz(zeta) * (EcF(rho) - EcP(rho));
           exc += rho * zk * grid_w_->pointer()[p];

        }else{
                double zk = 0.0;
                exc += 0.0;
             }
    }
    return exc;
}

// Conv. (I)
double MCPDFTSolver::EC_LYP_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                           std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB) {
    double tol = 1.0e-20;

    const double A = 0.04918;
    const double B = 0.132;
    const double C = 0.2533;
    const double Dd = 0.349;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_abp = SIGMA_AB->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double sigmaaa = sigma_aap[p];
        double sigmaab = sigma_abp[p];
        double sigmabb = sigma_bbp[p];
        double sigma = sigmaaa + sigmabb + 2.0 * sigmaab;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              double zk = 0.0;

           }else if ( rhob < tol ){

                    double zk = 0.0;
           }else{

                sigmabb = std::max(0.0,sigmabb);
                sigmaaa = std::max(0.0,sigmaaa);
                double rho_2    = pow(rho,2.0);
                double rho_a_2  = pow(rhoa,2.0);
                double rho_b_2  = pow(rhob,2.0);
                double rho_m13  = pow(rho,-1.0/3.0);
                double rho_a_83 = pow(rhoa,8.0/3.0);
                double rho_b_83 = pow(rhob,8.0/3.0);

                double f = pow( 1.0 + Dd * rho_m13, -1.0);

                double omega = exp(-C * rho_m13) * f * pow(rho,-11.0/3.0);

                double delta = (C + Dd * f) * rho_m13;

                double zk  = -4.0 * A  * rhoa * rhob * f / rho  - A * B * omega * (rhoa * rhob * (36.462398978764777098 * (rho_a_83 + rho_b_83)
                           + (2.6111111111111111111 - 0.38888888888888888889 * delta) * sigma - (2.5 - 0.055555555555555555556 * delta) * (sigmaaa + sigmabb)
                           - 0.11111111111111111111 * (delta - 11.0) * ((rhoa * sigmaaa + rhob * sigmabb) / rho ) ) - 0.66666666666666666667 * rho_2 * sigma
                           + (0.66666666666666666667 * rho_2 - rho_a_2) * sigmabb + (0.66666666666666666667 * rho_2 - rho_b_2) * sigmaaa);

                exc += zk * grid_w_->pointer()[p];
           }
           }else{
                double zk = 0.0;
                exc += zk * grid_w_->pointer()[p];
                }
    }
    return exc;
}

double MCPDFTSolver::EC_PBE(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                            std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB){

    double tol = 1.0e-20;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_abp = SIGMA_AB->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];

        double sigmaaa = sigma_aap[p];
        double sigmaab = sigma_abp[p];
        double sigmabb = sigma_bbp[p];

        double rho = rhoa + rhob;
        double sigma = sigmaaa + sigmabb + 2.0 * sigmaab;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;
              double t2 = 1.0 / rhob;
              double t3 = pow(t2 ,1.0/3.0);
              double t6 = pow(t2 ,1.0/6.0);
              double t9 = sqrt(t2);
              double t11 = t3 * t3;
              double t17 = log(1.0 + 0.3216395899738507e2 / (0.1112037486309468e2 * t6 + 0.3844746237447211e1 * t3 + 0.1644733775567609e1
                         * t9 + 0.2405871291288192 * t11));
              double t18 = (1.0 + 0.1274696188700087 * t3) * t17;
              double t20 = pow(rhob ,2.0);
              double t21 = pow(rhob ,1.0/3.0);
              double t23 = 1.0 / t21 / t20;
              double t26 = exp(0.2000000587336264e1 * t18);
              double t27 = t26 - 1.0;
              double t31 = 0.2162211495206379 / t27 * sigmabb * t23;
              double t33 = pow(t27 ,2.0);
              double t35 = pow(sigmabb ,2.0);
              double t37 = pow(t20 ,2.0);
              double t38 = pow(t21 ,2.0);
              double t49 = log(1.0 + 0.2162211495206379 * sigmabb * t23 * (1.0 + t31) /(1.0 + t31 + 0.4675158550002605e-1 / t33 * t35 / t38 / t37));
              double zk = rhob * (-0.310907e-1 * t18 + 0.1554534543482745e-1 * t49);

              exc += rho * zk * grid_w_->pointer()[p];

           }else if ( rhob < tol ){

                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;
                    double t2 = 1.0 / rhoa;
                    double t3 = pow(t2 ,1.0/3.0);
                    double t6 = pow(t2 ,1.0/6.0);
                    double t9 = sqrt(t2);
                    double t11 = t3 * t3;
                    double t17 = log(1.0 + 0.3216395899738507e2 / (0.1112037486309468e2 * t6 + 0.3844746237447211e1 * t3 + 0.1644733775567609e1
                               * t9 + 0.2405871291288192 * t11));
                    double t18 = (1.0 + 0.1274696188700087 * t3) * t17;
                    double t20 = pow(rhoa ,2.0);
                    double t21 = pow(rhoa ,1.0/3.0);
                    double t23 = 1.0 / t21 / t20;
                    double t26 = exp(0.2000000587336264e1 * t18);
                    double t27 = t26 - 1.0;
                    double t31 = 0.2162211495206379 / t27 * sigmaaa * t23;
                    double t33 = pow(t27 ,2.0);
                    double t35 = pow(sigmaaa ,2.0);
                    double t37 = pow(t20 ,2.0);
                    double t38 = pow(t21 ,2.0);
                    double t49 = log(1.0 + 0.2162211495206379 * sigmaaa * t23 * (1.0 + t31) /(1.0 + t31 + 0.4675158550002605e-1 / t33 * t35 / t38 / t37));
                    double zk = rhoa * (-0.310907e-1 * t18 + 0.1554534543482745e-1 * t49);

                    exc += rho * zk * grid_w_->pointer()[p];

           }else{

                double t4 = 1/rho;
                double t5 = pow( t4, 1.0/3.0);
                double t7 = 1.0 + 0.1325688999052018 * t5;
                double t8 = pow(t4, 1.0/6.0);
                double t11 = sqrt(t4);
                double t13 = pow(t5 ,2.0);
                double t15 = 0.598255043577108e1 * t8 + 0.2225569421150687e1 * t5 + 0.8004286349993634 * t11 + 0.1897004325747559 * t13;
                double t18 = 1.0 + 0.1608197949869254e2 / t15;
                double t19 = log(t18);
                double t21 = 0.621814e-1 * t7 * t19;
                double t23 = 1.0 + 0.6901399211255825e-1 * t5;
                double t28 = 0.8157414703487641e1 * t8 + 0.2247591863577616e1 * t5 + 0.4300972471276643 * t11 + 0.1911512595127338 * t13;
                double t31 = 1.0 + 0.2960874997779344e2 / t28;
                double t32 = log(t31);
                double t33 = t23 * t32;
                double t35 = rhoa - 1.0 * rhob;
                double t36 = t35 * t4;  // zeta
                double t37 = 1.0 + t36;
                double t38 = pow(t37 ,1.0/3.0);
                double t41 = 1.0 - t36;
                double t42 = pow(t41 ,1.0/3.0);
                double t44 = t38 * t37 + t42 * t41 - 2.0;
                double t45 = pow(t35 ,2.0);
                double t46 = pow(t45 ,2.0);
                double t47 = pow(rho ,2.0);
                double t48 = pow(t47 ,2.0);
                double t49 = 1.0 / t48;
                double t50 = t46 * t49;
                double t52 = 1.0 - t50;
                double t55 = 0.37995525e-1 * t33 * t44 * t52;
                double t57 = 1.0 + 0.1274696188700087 * t5;
                double t62 = 0.1112037486309468e2 * t8 + 0.3844746237447211e1 * t5 + 0.1644733775567609e1 * t11 + 0.2405871291288192 * t13;
                double t65 = 1.0 + 0.3216395899738507e2 / t62;
                double t66 = log(t65);
                double t69 = -0.310907e-1 * t57 * t66 + t21;
                double t70 = t69 * t44;
                double t72 = 0.1923661050931536e1 * t70 * t50;
                double t73 = pow(t38 ,2.0);
                double t75 = pow(t42 ,2.0);
                double t77 = 0.5 * t73 + 0.5 * t75;
                double t78 = pow(t77 ,2.0);
                double t79 = t78 * t77;
                double t80 = 1.0 / t78;
                double t81 = sigma * t80;
                double t82 = pow(rho ,1.0/3.0);
                double t84 = 1.0 / t82 / t47;
                double t85 = -t21 + t55 + t72;
                double t86 = 1.0 / t79;
                double t89 = exp(-0.3216396844291482e2 * t85 * t86);
                double t90 = t89 - 1.0;
                double t91 = 1.0 / t90;
                double t92 = t91 * sigma;
                double t93 = t80 * t84;
                double t95 = 0.1362107888567592 * t92 * t93;
                double t96 = 1.0 + t95;
                double t98 = pow(t90 ,2.0);
                double t99 = 1.0 / t98;
                double t100 = pow(sigma ,2.0);
                double t101 = t99 * t100;
                double t102 = pow(t78 ,2.0);
                double t103 = 1.0 / t102;
                double t104 = pow(t82 ,2.0);
                double t106 = 1.0 / t104 / t48;
                double t107 = t103 * t106;
                double t110 = 1.0 + t95 + 0.1855337900098064e-1 * t101 * t107;
                double t111 = 1.0 / t110;
                double t115 = 1.0 + 0.1362107888567592 * t81 * t84 * t96 * t111;
                double t116 = log(t115);
                double t118 = 0.310906908696549e-1 * t79 * t116;
                double zk = -t21 + t55 + t72 + t118;

                exc += rho * zk * grid_w_->pointer()[p];
           }
           }else{
                double zk = 0.0;
                exc += rho * zk * grid_w_->pointer()[p];
                }
    }
    return exc;
}

// From J.P. Perdew, and Y. Wang, Phys. Rev. B 45, 13244 (1992)
double MCPDFTSolver::EC_PW92_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B){

    const double pa = 1.0;
    const double Aa = 0.0168869;
    const double a1a = 0.11125;
    const double b1a = 10.357;
    const double b2a = 3.6231;
    const double b3a = 0.88026;
    const double b4a = 0.49671;
    const double pe = 1.0;
    const double c0p = 0.0310907;
    const double a1p = 0.21370;
    const double b1p = 7.5957;
    const double b2p = 3.5876;
    const double b3p = 1.6382;
    const double b4p = 0.49294;
    const double c0f = 0.01554535;
    const double a1f = 0.20548;
    const double b1f = 14.1189;
    const double b2f = 6.1977;
    const double b3f = 3.3662;
    const double b4f = 0.62517;
    const double d2Fz = 1.7099209341613656173;

    double tol = 1.0e-20;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    auto Fz = [](double ZETA) -> double {

              double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
              return dum;
    };

    auto G = [](double r, double T, double a1, double b1, double b2, double b3, double b4, double p) -> double {

             double dum = -2.0 * T * (1.0 + a1 * r) * log(1.0 + 0.5 * pow(T * (b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0)) ,-1.0));
             return dum;

    };

    auto Ac = [=](double r) -> double {

              double temp = -G(r,Aa,a1a,b1a,b2a,b3a,b4a,pa);
              return temp;
    };

    auto EcP = [=](double r) -> double {

               double dum = G(r,c0p,a1p,b1p,b2p,b3p,b4p,pe);
               return dum;
    };

    auto EcF = [=](double r) -> double {

               double dumm = G(r,c0f,a1f,b1f,b2f,b3f,b4f,pe);
               return dumm;
    };

    auto Ec = [=](double r, double ZETA) -> double {

              double dum = EcP(r) + ( Ac(r) * Fz(ZETA) * (1.0 - pow(ZETA ,4.0)) ) / d2Fz + ( EcF(r) - EcP(r) ) * Fz(ZETA) * pow(ZETA ,4.0);
              return dum;
    };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = (rhoa - rhob) / rho;
        double rs =  pow( 3.0 / ( 4.0 * M_PI * rho) , 1.0/3.0 );

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              zeta = 1.0;

           }else if ( rhob < tol ){

                    rho = rhoa;
                    zeta = 1.0;

           }
           double zk = Ec(rs,zeta);
           exc += rho * zk * grid_w_->pointer()[p];

        }else{
                double zk = 0.0;
                exc += rho * zk * grid_w_->pointer()[p];
             }
    }
    return exc;
}

// Conv (I)
double MCPDFTSolver::EC_PBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                           std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB) {
    const double pa = 1.0;
    const double Aa = 0.0168869;
    const double a1a = 0.11125;
    const double b1a = 10.357;
    const double b2a = 3.6231;
    const double b3a = 0.88026;
    const double b4a = 0.49671;
    const double pe = 1.0;
    const double c0p = 0.0310907;
    const double a1p = 0.21370;
    const double b1p = 7.5957;
    const double b2p = 3.5876;
    const double b3p = 1.6382;
    const double b4p = 0.49294;
    const double c0f = 0.01554535;
    const double a1f = 0.20548;
    const double b1f = 14.1189;
    const double b2f = 6.1977;
    const double b3f = 3.3662;
    const double b4f = 0.62517;
    const double d2Fz = 1.7099209341613656173;
    const double BETA = 0.06672455060314922;
    const double GAMMA = 0.0310906908696549;

    double tol = 1.0e-20;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_abp = SIGMA_AB->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

        auto Fi = [=](double ZETA) -> double {

                  double dumm = 0.5 * (pow((1.0 + ZETA) ,2.0/3.0 ) + pow((1.0 - ZETA) ,2.0/3.0));
                  return dumm;
        };

        auto kF = [](double RHO) -> double {

                  double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
                  return dum;
        };

        auto ks = [=](double RHO) -> double {

                  double temp = sqrt(4.0 * kF(RHO) / M_PI);
                  return temp;
        };

        auto Fz = [](double ZETA) -> double {

                  double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
                  return dum;
        };

        auto t = [=](double RHO, double SIGMA, double ZETA) -> double {

                 double temp = sqrt(SIGMA) / (2.0 * ks(RHO) * Fi(ZETA) * RHO);
                 return temp;
        };

        auto G = [](double r, double T, double a1, double b1, double b2, double b3, double b4, double p) -> double {

                 double dum = -2.0 * T * (1.0 + a1 * r) * log(1.0 + 0.5 * pow(T * (b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0)) ,-1.0));
                 return dum;

        };

        auto Ac = [=](double r) -> double {

                  double temp = -G(r,Aa,a1a,b1a,b2a,b3a,b4a,pa);
                  return temp;
        };

        auto EcP = [=](double r) -> double {

                   double dum = G(r,c0p,a1p,b1p,b2p,b3p,b4p,pe);
                   return dum;
        };

        auto EcF = [=](double r) -> double {

                   double dumm = G(r,c0f,a1f,b1f,b2f,b3f,b4f,pe);
                   return dumm;
        };

        auto Ec = [=](double r, double ZETA) -> double {

                  double dum = EcP(r) + ( Ac(r) * Fz(ZETA) * (1.0 - pow(ZETA ,4.0)) ) / d2Fz + ( EcF(r) - EcP(r) ) * Fz(ZETA) * pow(ZETA ,4.0);
                  return dum;
        };

        auto A = [=](double r, double ZETA) -> double {

                 double dum = (BETA/GAMMA) * pow( exp(-Ec(r,ZETA) / (pow(Fi(ZETA),3.0) * GAMMA)) - 1.0, -1.0);
                 return dum;
        };

        auto H = [=](double RHO, double SIGMA, double r, double ZETA) -> double {

                 double temp = pow(Fi(ZETA),3.0) * GAMMA * log(1.0 + (BETA/GAMMA) * pow(t(RHO,SIGMA,ZETA) ,2.0) * (1.0 + A(r,ZETA) * pow(t(RHO,SIGMA,ZETA),2.0))
                             / (1.0 + A(r,ZETA) * pow(t(RHO,SIGMA,ZETA),2.0) + pow(A(r,ZETA),2.0) * pow(t(RHO,SIGMA,ZETA),4.0)));
                 return temp;
        };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = (rhoa - rhob) / rho;
        double rs =  pow( 3.0 / ( 4.0 * M_PI * rho) , 1.0/3.0 );
        double sigmaaa = sigma_aap[p];
        double sigmaab = sigma_abp[p];
        double sigmabb = sigma_bbp[p];
        double sigma = sigmaaa + sigmabb + 2.0 * sigmaab;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;
              zeta = 1.0;

           }else if ( rhob < tol ){

                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;
                    zeta = 1.0;

           }else/* if (!(rhoa < tol) && !(rhob < tol) ) */{
                   
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigmabb = std::max(0.0,sigmabb);
                    sigma = sigmaaa + sigmabb + 2.0 * sigmaab;

           }
           double zk = H(rho,sigma,rs,zeta) + Ec(rs,zeta);
           exc += rho * zk * grid_w_->pointer()[p];
          
        }else{
                double zk = 0.0;
                exc += rho * zk * grid_w_->pointer()[p];
             }
    }
    return exc;
}

}} // End namespaces
