/*
 * @BEGIN LICENSE
 *
 * mcpdft by Psi4 Developer, a plugin to:
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

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "psi4/libpsi4util/libpsi4util.h"

#include "psi4/libqt/qt.h"

// jk object
#include "psi4/libfock/jk.h"

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

// for potential object
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

// mcpdft 
#include "mcpdft_solver.h"

// blas
#include "blas.h"
#include "blas_mangle.h"

// for reading 2RDM
#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>

#include <psi4/libpsi4util/PsiOutStream.h>

using namespace psi;
using namespace fnocc;

namespace psi{ namespace mcpdft {

// the MCPDFTSolver class derives from the Wavefunction class and inherits its members
MCPDFTSolver::MCPDFTSolver(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options_):
    Wavefunction(options_){

    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

MCPDFTSolver::~MCPDFTSolver() {

}

void MCPDFTSolver::print_banner() {
   outfile->Printf("\n\n");
   outfile->Printf( "        ********************************************************************\n");
   outfile->Printf( "        *                                                                  *\n");
   outfile->Printf( "        *    MCPDFT: Multiconfigurational Pair-Density Functional Theory   *\n");
   outfile->Printf( "        *                                                                  *\n");
   outfile->Printf( "        *    Mohammad Mostafanejad and A. Eugene DePrince III              *\n");
   outfile->Printf( "        *                                                                  *\n");
   outfile->Printf( "        ********************************************************************\n");
   outfile->Printf("\n\n");

   outfile->Printf("\n\n");
   outfile->Printf("        The following paper should be cited when using MCPDFT:\n");
   outfile->Printf("\n");
   outfile->Printf("        M. Mostafanejad, and A. E. DePrince III,\n");
   outfile->Printf("        J. Chem. Theory Comput. 15, 290-302 (2019).\n");
   outfile->Printf("\n");
   outfile->Printf("        URL: https://doi.org/10.1021/acs.jctc.8b00988\n");
   outfile->Printf("\n");
   outfile->Printf("\n\n");
}

// initialize members of the MCPDFTSolver class
void MCPDFTSolver::common_init() {

    is_low_memory_ = false;

    reference_energy_ = Process::environment.globals["V2RDM TOTAL ENERGY"];
    
    if (  options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT"
       || options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT" ) {
       
       // calculating mp2 energy for double-hybrids
       mp2_corr_energy_ = Process::environment.globals["MP2 CORRELATION ENERGY"];
    }
    
    shallow_copy(reference_wavefunction_);

    // number of alpha electrons
    nalpha_   = reference_wavefunction_->nalpha();

    // number of beta electrons
    nbeta_    = reference_wavefunction_->nbeta();

    // number of alpha electrons per irrep
    nalphapi_ = reference_wavefunction_->nalphapi();

    // number of beta electrons per irrep
    nbetapi_  = reference_wavefunction_->nbetapi();

    // number of doubly occupied orbitals per irrep
    doccpi_   = reference_wavefunction_->doccpi();

    // number of singly occupied orbitals per irrep
    soccpi_   = reference_wavefunction_->soccpi();

    // number of frozen core orbitals per irrep
    frzcpi_   = reference_wavefunction_->frzcpi();

    // number of frozen virtual orbitals per irrep
    frzvpi_   = reference_wavefunction_->frzvpi();

    // number of molecular orbials per irrep
    nmopi_    = reference_wavefunction_->nmopi();

    // number of symmetry orbials per irrep
    nsopi_    = reference_wavefunction_->nsopi();

    // number of irreducible representations
    nirrep_   = reference_wavefunction_->nirrep();

    // total number of symmetry orbitals
    nso_      = reference_wavefunction_->nso();

    // total number of molecular orbitals
    nmo_      = reference_wavefunction_->nmo();

    // if (  options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT"
    //    || options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT" ) {

    //    int nrstdocc = 0;
    //    int nelec = 0;
    //    for (int h = 0; h < nirrep_; h++){
    //        nrstdocc += doccpi_[h];
    //        nelec += nalphapi_[h];
    //        nelec += nbetapi_[h];
    //        printf("docc nelec = %i %i\n",nrstdocc,nelec);
    //    }
    //    if (nelec != nrstdocc) throw PsiException("All electrons should be frozen for double-hybrids!\n \
    //         The reference wave function should be single-determinant.\n",__FILE__,__LINE__);
    // }

    // grab the molecule from the reference wave function
    molecule_ = reference_wavefunction_->molecule();

    // SO/MO transformation matrices
    Ca_ = std::shared_ptr<Matrix>(reference_wavefunction_->Ca());
    Cb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Cb());

    // overlap integrals
    S_  = std::shared_ptr<Matrix>(reference_wavefunction_->S());

    // SO-basis Fock matrices
    Fa_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fa());
    Fb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fb());

    // SO-basis density matrices
    Da_ = std::shared_ptr<Matrix>(reference_wavefunction_->Da());
    Db_ = std::shared_ptr<Matrix>(reference_wavefunction_->Db());

    // orbital energies
    epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // set the wavefunction name
    name_ = "DFT";

    // restricted orbitals, unrestricted rdms
    same_a_b_orbs_ = true;
    same_a_b_dens_ = false;

    // symmetry of orbitals:
    symmetry_ = (int*)malloc(nmo_*sizeof(int));
    memset((void*)symmetry_,'\0',nmo_*sizeof(int));
    int count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nsopi_[h]; i++) {
            symmetry_[count++] = h;
        }
    }

    // pitzer offset
    pitzer_offset_ = (int*)malloc(nirrep_*sizeof(int));
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset_[h] = count;
        count += nsopi_[h];
    }

    // evaluate basis function values on a grid:

    scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    potential_ = scfwfn->V_potential();

    // phi matrix (real-space <- AO basis mapping)
    // since grid is stored in blocks, we need to build a full phi matrix
    // from the individual blocks:

    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");
    int nao = primary->nbf();
    std::shared_ptr<Matrix> Da_ao (new Matrix(nao,nao));
    std::shared_ptr<Matrix> Db_ao (new Matrix(nao,nao));

    points_func_ = potential_->properties()[0];
    points_func_->set_pointers(Da_ao,Db_ao);

    // determine number of grid points
    int nblocks = potential_->nblocks();
    phi_points_    = 0;
    max_functions_ = 0;
    max_points_    = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        phi_points_ += npoints;
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();
        if ( nlocal > max_functions_ ) max_functions_ = nlocal;
        if ( npoints > max_points_ )   max_points_    = npoints;
    }

    // what is the derivative level?
    deriv_   = functional->deriv();
    is_gga_  = functional->is_gga();
    is_meta_ = functional->is_meta();
    is_unpolarized_ = functional->is_unpolarized();

    std::vector < std::shared_ptr<Matrix> > JK = BuildJK();

    double caa = Da_->vector_dot(JK[0]);
    double cab = Da_->vector_dot(JK[1]);
    double cba = Db_->vector_dot(JK[0]);
    double cbb = Db_->vector_dot(JK[1]);

    coulomb_energy_ = 0.5 * ( caa + cab + cba + cbb );

    // hf_ex_energy_ = 0.0;
    // lr_ex_energy_ = 0.0;
    if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {

       // HF exchange energy should be computed using K object
       double kaa = Da_->vector_dot(JK[2]);
       double kbb = Db_->vector_dot(JK[3]);
 
       hf_ex_energy_ = -0.5 * (kaa + kbb);
    
       if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
          && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {

          // long range (LR) exchange energy calculated using JK object
          double wkaa = Da_->vector_dot(JK[4]);
          double wkbb = Db_->vector_dot(JK[5]);
 
          lr_ex_energy_ = -0.5 * (wkaa + wkbb);
       }
    }
    
    /* ==============================================
       estimating the memory requirement for MCPDFT
       ============================================== */

    outfile->Printf("\n"); 
    outfile->Printf("    ==> Memory requirements <==\n");
    outfile->Printf("\n");
    
    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // assume 500 MB is already allocated by psi4
    long int n_extra       = 500 * 1024 * 1024 / 8;
    long int n_mem         = n_extra;
    long int available_mem = 0;
    long int n_init        = 0;
    long int n_phiAO       = 0;
    long int n_grids       = 0;
    long int n_transf      = 0;
    long int n_rho         = 0;
    long int n_pi          = 0;
    long int n_R           = 0;
    long int n_transl      = 0;
    long int n_df          = 0;

    // default basic vectors and matrices initialized using common_init()
    n_init  = nso_ * nso_;                                 // S_ matrix
    n_init += 2 * nso_ * nso_;                             // Ca_ and Cb_ matrices
    n_init += 2 * nso_ * nso_;                             // Fa_ and Fb_ matrices
    n_init += 2 * nmo_;                                    // epsilon_a_ and epsilon_b_ vectors
    n_mem += n_init;

    // phiAO, phiAO_x, phiAO_y, phiAO_z
    n_phiAO = phi_points_ * nso_;
    if ( is_gga_ || is_meta_ ) { 
       n_phiAO   += 3 * phi_points_ * nso_;
    }
    n_mem += n_phiAO;

    // required memory for grids x, y, and z and weights w vectors
    n_grids = 4 * phi_points_;
    n_mem += n_grids;

    // memory needed for AO->MO transformation
    n_transf = nirrep_ * phi_points_;                     // phi_points_list vector
    n_transf += phi_points_ * nso_;                       // super_phi_ matrix
    n_transf += nso_ * nso_;                              // AO2SO matrix in TransformPhiMatrixAOMO()
    n_transf += phi_points_ * nso_;                       // temporary matrix called three times in TransformPhiMatrixAOMO()
    if ( is_gga_ || is_meta_ ) { 
        n_transf += 3 * phi_points_ * nso_;               // for super_phi_x_, _y_ and _z_ gradient matrices
    }
    n_mem += n_transf;

    // memory needed for storing rho and rho' matrices
    n_rho  = 3 * phi_points_;                             // rho_a_, rho_b_ and rho_ vectors
    if ( is_gga_ || is_meta_ ) { 
        n_rho += 3 * phi_points_;                         // rho_a_x_, rho_a_y_ and rho_a_z_ gradient vectors
        n_rho += 3 * phi_points_;                         // rho_b_x_, rho_b_y_ and rho_b_z_ gradient vectors
        n_rho += 3 * phi_points_;                         // sigma_aa_, sigma_ab_ and sigma_bb_ vectors
    }
    n_mem += n_rho;

    // memory needed for storing pi vector
    n_pi   = phi_points_;                                 // pi_ vector
    if ( is_gga_ || is_meta_ ) { 
        n_pi  += 3 * phi_points_;                         // pi_x_, pi_y_ and pi_z_ gradient vectors
    }
    n_mem += n_pi;
 
    // memory needed for storing on-top ratio vector
    n_R    = phi_points_;
    n_mem += n_R;
    
    // memory needed for (full-) translation step
    n_transl  = 2 * phi_points_;                          // (f)tr_rho_a_, (f)tr_rho_b_ vectors
    if ( is_gga_ || is_meta_ ) { 
        n_transl += 3 * phi_points_;                      // (f)tr_rho_a_x_, (f)tr_rho_a_y_ and (f)tr_rho_a_z_ gradient vectors
        n_transl += 3 * phi_points_;                      // (f)tr_rho_b_x_, (f)tr_rho_b_y_ and (f)tr_rho_b_z_ gradient vectors
        n_transl += 3 * phi_points_;                      // (f)tr_sigma_aa_, (f)tr_sigma_ab_ and (f)tr_sigma_bb_ vectors
    }
    n_mem += n_transl;

    if ( n_mem * 8.0 > (double)memory_ ) {
        outfile->Printf("\n"); 
        outfile->Printf("      <<< Warning >>> \n");
        outfile->Printf("\n");
        outfile->Printf("      Not enough memory. Switching to low-memory algorithm.\n");
        outfile->Printf("      For optimal performance, increase the available memory\n");
        outfile->Printf("      by at least %7.2lf mb.\n",(8.0 * n_mem - memory_)/1024.0/1024.0);
        outfile->Printf("\n");

        is_low_memory_ = true;
        n_mem -= n_transf;
        n_mem -= n_phiAO;
        if ( n_mem * 8.0 > (double)memory_ ) {
            outfile->Printf("\n");
            outfile->Printf("      Wow.  On second thought, there isn't even enough memory for the low-memory algorithm!\n");
            outfile->Printf("\n");
            throw PsiException("Not enough memory.",__FILE__,__LINE__); fflush(stdout);
        }
    }

    outfile->Printf("    =========================================================\n");
    outfile->Printf("    memory specified by the user:                %7.2lf MB\n",(double)memory_ / 1024.0 / 1024.0);
    outfile->Printf("    ---------------------------------------------------------\n");
    outfile->Printf("    initialization:                              %7.2lf MB\n",n_init  * sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    grid-points and weights:                     %7.2lf MB\n",n_grids * sizeof(double)   / 1024.0 / 1024.0);
    if ( !is_low_memory_ ) {
        outfile->Printf("    phi & phi' (AO):                             %7.2lf MB\n",n_phiAO * sizeof(double)   / 1024.0 / 1024.0);
        outfile->Printf("    AO->SO transformation:                       %7.2lf MB\n",n_transf* sizeof(double)   / 1024.0 / 1024.0);
    }
    outfile->Printf("    rho and rho':                                %7.2lf MB\n",n_rho *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    on-top pair density:                         %7.2lf MB\n",n_pi  *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    on-top ratio:                                %7.2lf MB\n",n_R   *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    translation step:                            %7.2lf MB\n",n_transl * sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    extra:                                       %7.2lf MB\n",n_extra *  sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    ---------------------------------------------------------\n");
    outfile->Printf("    total memory for MCPDFT:                     %7.2lf MB\n",n_mem *    sizeof(double)   / 1024.0 / 1024.0);
    outfile->Printf("    =========================================================\n");

    // memory available after allocating all we need for MCPDFT (deleting the temporary matrices and vectors)
    n_mem -= nirrep_ * phi_points_;                // phi_points_list vector
    n_mem -= nso_ * nso_;                          // AO2SO matrix in TransformPhiMatrixAOMO()
    if ( !is_low_memory_ ) {
        n_mem -= n_phiAO;                          // ao-basis phi matrices
        n_mem -= phi_points_ * nso_;               // temporary matrix called three times in TransformPhiMatrixAOMO()
    }

    available_memory_ = memory_ - n_mem * 8L;

    outfile->Printf("\n");
    outfile->Printf("    available memeory for building JK object:    %7.2lf MB\n",(double)available_memory_ / 1024.0 / 1024.0);
    outfile->Printf("\n");

    grid_x_      = std::shared_ptr<Vector>(new Vector("GRID X",phi_points_));
    grid_y_      = std::shared_ptr<Vector>(new Vector("GRID Y",phi_points_));
    grid_z_      = std::shared_ptr<Vector>(new Vector("GRID Z",phi_points_));
    grid_w_      = std::shared_ptr<Vector>(new Vector("GRID W",phi_points_));

    if ( is_low_memory_ ) {

        GetGridInfo();

    }else {
        outfile->Printf("\n"); 
        outfile->Printf("    ==> Build Phi and Phi' matrices ..."); fflush(stdout);

        std::shared_ptr<Matrix> super_phi_ao (new Matrix("SUPER PHI (AO)",phi_points_,nso_));

        std::shared_ptr<Matrix> super_phi_x_ao;
        std::shared_ptr<Matrix> super_phi_y_ao;
        std::shared_ptr<Matrix> super_phi_z_ao;

        if ( is_gga_ || is_meta_ ) {

            super_phi_x_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X (AO)",phi_points_,nso_));
            super_phi_y_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y (AO)",phi_points_,nso_));
            super_phi_z_ao = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z (AO)",phi_points_,nso_));
            
        }

        // build phi matrix and derivative phi matrices (ao basis)
        BuildPhiMatrixAO("PHI",super_phi_ao);

        if ( is_gga_ || is_meta_ ) {

            BuildPhiMatrixAO("PHI_X",super_phi_x_ao);
            BuildPhiMatrixAO("PHI_Y",super_phi_y_ao);
            BuildPhiMatrixAO("PHI_Z",super_phi_z_ao);

        }

        // transform the orbital label in phi matrix and derivative phi matrix to the MO basis

        Dimension phi_points_list = Dimension(nirrep_,"phi_points");
        for (int h = 0; h < nirrep_; h++) {
            phi_points_list[h] = phi_points_;
        }

        super_phi_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI",phi_points_list,nsopi_));

        TransformPhiMatrixAOMO(super_phi_ao,super_phi_);

        if ( is_gga_ || is_meta_ ) {

            super_phi_x_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X",phi_points_list,nsopi_));
            super_phi_y_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y",phi_points_list,nsopi_));
            super_phi_z_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z",phi_points_list,nsopi_));

            TransformPhiMatrixAOMO(super_phi_x_ao,super_phi_x_);
            TransformPhiMatrixAOMO(super_phi_y_ao,super_phi_y_);
            TransformPhiMatrixAOMO(super_phi_z_ao,super_phi_z_);

        }
        outfile->Printf(" Done. <==\n"); 
    }

}

void MCPDFTSolver::GetGridInfo() {

    int nblocks = potential_->nblocks();

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func_->basis_value("PHI")->pointer();

        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();

        for (int p = 0; p < npoints; p++) {

            grid_x_->pointer()[phi_points_ + p] = x[p];
            grid_y_->pointer()[phi_points_ + p] = y[p];
            grid_z_->pointer()[phi_points_ + p] = z[p];
            grid_w_->pointer()[phi_points_ + p] = w[p];

        }
        phi_points_ += npoints;
    }
}

void MCPDFTSolver::BuildPhiMatrixAO(std::string phi_type, std::shared_ptr<Matrix> myphi) {

    int nblocks = potential_->nblocks();

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func_->basis_value(phi_type)->pointer();

        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();

        // Doing some test to see everything including Pi etc is correct on 
        // Molcas' grid points through comparison.

        // std::ifstream dataIn;
       
        // dataIn.open("H2.grids_test");
        // 
        // if (!dataIn)
        //    std::cout << "Error opening file.\n";
        // else { 
        //      int p = 0;        
        //      while (!dataIn.eof()){
        //    
        //            dataIn >> x[p];
        //            dataIn >> y[p];    
        //            dataIn >> z[p];
        //            p++;
        //      }        
        // }
        // dataIn.close(); 

        // for (int p = 0; p < npoints; p++) {

        //     outfile->Printf("\n     y[");
        //     outfile->Printf("%d",p);
        //     outfile->Printf("] = %20.7lf\n",y[p]);
        // }

        for (int p = 0; p < npoints; p++) {

            grid_x_->pointer()[phi_points_ + p] = x[p];
            grid_y_->pointer()[phi_points_ + p] = y[p];
            grid_z_->pointer()[phi_points_ + p] = z[p];
            grid_w_->pointer()[phi_points_ + p] = w[p];

            for (int nu = 0; nu < nlocal; nu++) {
                int nug = function_map[nu];
                myphi->pointer()[phi_points_ + p][nug] = phi[p][nu];
            }
        }
        phi_points_ += npoints;
    }
}

void MCPDFTSolver::TransformPhiMatrixAOMO(std::shared_ptr<Matrix> phi_in, std::shared_ptr<Matrix> phi_out) {
        
    std::shared_ptr<Matrix> phi_temp (new Matrix(phi_out));

    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] != 0 )
           F_DGEMM('n','n',nsopi_[h],phi_points_,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],phi_in->pointer()[0],nso_,0.0,phi_temp->pointer(h)[0],nsopi_[h]);
    }

    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] != 0 )
           F_DGEMM('n','n',nsopi_[h],phi_points_,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],phi_temp->pointer(h)[0],nsopi_[h],0.0,phi_out->pointer(h)[0],nsopi_[h]);
    }
}

double MCPDFTSolver::compute_energy() {
    
    /* ========================================================================================
       Calculation of the range-separated global (double-)hybrid formalis is based
       on Eq. (12) of the reference: C. Kalai and J, Toulouse J. Chem. Phys. 148, 164105 (2018)
       ======================================================================================== */

    // read 1- and 2-RDM from disk and build rho(r), rho'(r), pi(r), and pi'(r)

    if ( options_.get_str("MCPDFT_REFERENCE") == "V2RDM" ) {

        ReadOPDM();
        ReadTPDM();

    }else if ( options_.get_str("MCPDFT_REFERENCE") == "CI" ) {

        if ( nirrep_ > 1 ) {
            throw PsiException("error, MCPDFT_REFERENCE = CI only works with symmetry c1",__FILE__,__LINE__);
        }

        // read 1-RDM
        ReadCIOPDM(Da_,"opdm_a.txt");
        ReadCIOPDM(Db_,"opdm_b.txt");

        // allocate memory 2-RDM
        double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
        memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

        // read 2-RDM
        ReadCITPDM(D2ab,"tpdm_ab.txt");

        // build alpha- and beta-spin densities and gradients (already built for MCPDFT_REFERENCE = V2RDM)
        outfile->Printf("\n"); 
        outfile->Printf("    ==> Build Rho <== \n ");

        BuildRho();

        // build on-top pair density (already built for MCPDFT_REFERENCE = V2RDM)
        outfile->Printf("\n");
        outfile->Printf("    ==> Build Pi <==\n\n");

        BuildPi(D2ab);
        free(D2ab);

    }else {
        throw PsiException("invalid MCPDFT_REFERENCE type",__FILE__,__LINE__);
    }

    // build R(r) = 4 * Pi(r) / rho(r)
    outfile->Printf("    ==> Build the on-top ratio R ...");
    Build_R();
    outfile->Printf(" Done. <==\n\n");

    // translate the alpha and beta densities and their corresponding gradients
    if ( options_.get_str("MCPDFT_TRANSLATION_TYPE") == "REGULAR") {

        outfile->Printf("    ==> Regular translation of densities and/or density gradients ...\n");
        Translate();    
        outfile->Printf("    ... Done. <==\n\n");

    }else {

        outfile->Printf("    ==> Full translation of densities and/or density gradients ...\n");
        Fully_Translate();    
        outfile->Printf("    ... Done. <==\n\n");

    }

    // calculate the on-top energy

    outfile->Printf("    ==> Evaluate the on-top energy contribution <==\n");
    outfile->Printf("\n");

    /* ===================================================================================================
       calculate the complement short-range MCPDFT XC functional energy:
       E = min(Psi->N) { <Psi| T + Wee_LR(w) + lambda * Wee_SR(w) + Vne |Psi> + E_HXC_(w,lambda)[rho,pi] }
       =================================================================================================== */

    // initialization of the hybrid (lambda) and range-separation (omega) parameters
    double mcpdft_lambda    = 0.0;
    double mcpdft_omega     = 0.0;
    if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {
       
       if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
          && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {

          // extracting the omega value from input file (default w=0.0)
          mcpdft_omega = options_.get_double("MCPDFT_OMEGA");

          outfile->Printf("    ==> Range-separation parameter (omega) = %5.2lf ", mcpdft_omega);
          outfile->Printf("<==\n");

          // calculating the LR and SR two-electron energies
          lr_Vee_energy_ = RangeSeparatedTEE("LR");
          sr_Vee_energy_ = RangeSeparatedTEE("SR");

       }
       // extracting the lambda value from input file (default lambda=0.0)
       mcpdft_lambda = options_.get_double("MCPDFT_LAMBDA");
       
       outfile->Printf("    ==> Coupling parameter for hybrid MCPDFT (lambda) = %5.2lf ", mcpdft_lambda);
       outfile->Printf("<==\n");
    }

    // computing the exchange and correlation density functional energies
    double mcpdft_ex = 0.0;
    double mcpdft_ec = 0.0;
    if (  (options_.get_str("MCPDFT_METHOD") == "RS_MCPDFT") 
       || (options_.get_str("MCPDFT_METHOD") == "RS1H_MCPDFT")
       || (options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT") ) {

       if ( options_.get_str("MCPDFT_FUNCTIONAL") == "WPBE" ) {
 
          mcpdft_ex = EX_wPBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
          mcpdft_ec = EC_PBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_ab_, tr_sigma_bb_);

       }else if ( options_.get_str("MCPDFT_FUNCTIONAL") == "WB97X" ) {
             throw PsiException("The blame is on Marcus who has to implement this functional!",__FILE__,__LINE__);
       }else {
             throw PsiException("Sorry! The requested range-separated functional has not yet been implemented!",__FILE__,__LINE__);
       }
    }else if (  (options_.get_str("MCPDFT_METHOD") == "MCPDFT") 
             || (options_.get_str("MCPDFT_METHOD") == "1H_MCPDFT")
             || (options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT") ) {

             if ( options_.get_str("MCPDFT_FUNCTIONAL") == "SVWN" ) {

                mcpdft_ex = EX_LSDA(tr_rho_a_, tr_rho_b_);
                mcpdft_ec = EC_VWN3_RPA_III(tr_rho_a_, tr_rho_b_);

             }else if ( options_.get_str("MCPDFT_FUNCTIONAL") == "BLYP" ) {
                
                      mcpdft_ex = EX_B88_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
                      mcpdft_ec = EC_LYP_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_ab_, tr_sigma_bb_);
             
             }else if (  options_.get_str("MCPDFT_FUNCTIONAL") == "PBE" 
                      || options_.get_str("MCPDFT_FUNCTIONAL") == "REVPBE" ) {
 
                      mcpdft_ex = EX_PBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
                      mcpdft_ec = EC_PBE_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_ab_, tr_sigma_bb_);

             }else if ( options_.get_str("MCPDFT_FUNCTIONAL") == "BOP" ) {
                
                      mcpdft_ex = EX_B88_I(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
                      mcpdft_ec = EC_B88_OP(tr_rho_a_, tr_rho_b_, tr_sigma_aa_, tr_sigma_bb_);
             }else{
                  throw PsiException("Please choose a proper functional for MCDPFT",__FILE__,__LINE__);
             }
    }else {
          throw PsiException("Please choose a valid method for MCDPFT",__FILE__,__LINE__);
    }
    
    // evaluate the kinetic, potential, and coulomb energies
    outfile->Printf("    ==> Evaluate kinetic, potential, and coulomb energies <==\n");
    outfile->Printf("\n");

    // one-electron terms:
    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));

    // SharedMatrix ha (new Matrix(mints->so_potential()));
    // ha->add(mints->so_kinetic());
    // ha->transform(Ca_);

    // SharedMatrix hb (new Matrix(mints->so_potential()));
    // hb->add(mints->so_kinetic());
    // hb->transform(Cb_);

    // double one_electron_energy = Da_->vector_dot(ha) 
    //                            + Db_->vector_dot(hb);

    double one_electron_energy = 0.0;

    // kinetic-energy 
    SharedMatrix Ta (new Matrix(mints->so_kinetic()));
    Ta->transform(Ca_);
    
    SharedMatrix Tb (new Matrix(mints->so_kinetic()));
    Tb->transform(Cb_);
    
    double kinetic_energy = Da_->vector_dot(Ta) 
                          + Db_->vector_dot(Tb);

    // nuclear-attraction potential energy
    SharedMatrix Va (new Matrix(mints->so_potential()));
    Va->transform(Ca_);
    
    SharedMatrix Vb (new Matrix(mints->so_potential()));
    Vb->transform(Cb_);
    
    double nuclear_attraction_energy = Da_->vector_dot(Va) 
                                     + Db_->vector_dot(Vb);
    
    one_electron_energy  = nuclear_attraction_energy + kinetic_energy;

//printf("yay!\n\n"); fflush(stdout);
//    // Coulomb energy computed using J object
//    // std::vector < std::shared_ptr<Matrix> > JK = BuildJK();
//
//printf("yay!\n\n"); fflush(stdout);
//    double caa = Da_->vector_dot(JK[0]);
//printf("yay!\n\n"); fflush(stdout);
//    double cab = Da_->vector_dot(JK[1]);
//printf("yay!\n\n"); fflush(stdout);
//    double cba = Db_->vector_dot(JK[0]);
//printf("yay!\n\n"); fflush(stdout);
//    double cbb = Db_->vector_dot(JK[1]);
//printf("yay!\n\n"); fflush(stdout);
//
//    double coulomb_energy = 0.5 * ( caa + cab + cba + cbb );

    // classical nuclear repulsion energy
    double nuclear_repulsion_energy = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // two-electron energy from reference wavefunction:  < Psi|  r12^-1 | Psi >
    two_electron_energy_ = reference_energy_ - nuclear_repulsion_energy - one_electron_energy;


//     double hf_ex_energy = 0.0;
//     double lr_ex_energy = 0.0;
//     if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {
// 
//        // HF exchange energy should be computed using K object
//        double kaa = Da_->vector_dot(JK[2]);
//        double kbb = Db_->vector_dot(JK[3]);
//  
//        hf_ex_energy = -0.5 * (kaa + kbb);
//     
//        if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
//           && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {
// 
//           // long range (LR) exchange energy calculated using JK object
//           double wkaa = Da_->vector_dot(JK[4]);
//           double wkbb = Db_->vector_dot(JK[5]);
//  
//           lr_ex_energy = -0.5 * (wkaa + wkbb);
//        }
//     }
    // printf("Coulomb_energy    %20.12lf \n", coulomb_energy);
    // printf("hartree_ex_energy %20.12lf \n", hf_ex_energy);
    // printf("lr_ex_energy      %20.12lf \n", lr_ex_energy);
    // printf("J-K               %20.12lf \n",coulomb_energy+ 0.25* hf_ex_energy);
    // printf("J-wK              %20.12lf \n",coulomb_energy+ 0.75* lr_ex_energy);
    // printf("J-K-wk            %20.12lf \n",coulomb_energy+ 0.25* hf_ex_energy + 0.75* lr_ex_energy);
    // printf("two electron energies (ref, total, sr, lr) = %20.12lf %20.12lf %20.12lf %20.12lf \n\n",two_electron_energy_,sr_Vee_energy_+lr_Vee_energy_,sr_Vee_energy_,lr_Vee_energy_);

    // print total energy and its components
    outfile->Printf("    ==> Energetics <==\n");
    outfile->Printf("\n");

    // outfile->Printf("        nuclear repulsion energy =    %20.12lf\n",molecule_->nuclear_repulsion_energy());
    outfile->Printf("        nuclear repulsion energy  =         %20.12lf\n",nuclear_repulsion_energy);
    outfile->Printf("        nuclear attraction energy =         %20.12lf\n",nuclear_attraction_energy);
    outfile->Printf("        kinetic energy            =         %20.12lf\n",kinetic_energy);
    outfile->Printf("        one-electron energy       =         %20.12lf\n",one_electron_energy);
    if ( (options_.get_str("MCPDFT_METHOD") == "1H_MCPDFT") || (options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT") ) {
       outfile->Printf("        two-electron energy       =         %20.12lf\n",two_electron_energy_);
       outfile->Printf("        HF-exchange energy        =         %20.12lf\n",hf_ex_energy_);
       if ( options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT") {
          outfile->Printf("        MP2-correlation energy    =         %20.12lf\n",mp2_corr_energy_);
       }
    }else if( (options_.get_str("MCPDFT_METHOD") == "RS1H_MCPDFT") || (options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT") ) {
       outfile->Printf("        two-electron energy       =         %20.12lf\n",two_electron_energy_);
       outfile->Printf("        HF-exchange energy        =         %20.12lf\n",hf_ex_energy_);
       outfile->Printf("        LR-exchange enrgy        =         %20.12lf\n",lr_ex_energy_);
       if ( options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT") {
          outfile->Printf("        MP2-correlation energy    =         %20.12lf\n",mp2_corr_energy_);
       }
    }
    outfile->Printf("        classical coulomb energy  =         %20.12lf\n",coulomb_energy_);
    if ( options_.get_str("MCPDFT_REFERENCE") == "V2RDM") {
       outfile->Printf("        v2RDM-CASSCF energy contribution =  %20.12lf\n",
                      nuclear_repulsion_energy + one_electron_energy + coulomb_energy_);
    }else{
         outfile->Printf("        CASSCF energy contribution =        %20.12lf\n",
                        nuclear_repulsion_energy + one_electron_energy + coulomb_energy_);
    }
    outfile->Printf("        Ex                        =         %20.12lf\n",mcpdft_ex);
    outfile->Printf("        Ec                        =         %20.12lf\n",mcpdft_ec);
    outfile->Printf("        On-top energy =                     %20.12lf\n",mcpdft_ex + mcpdft_ec);
    outfile->Printf("\n");

    // double total_energy = nuclear_repulsion_energy + one_electron_energy + lmbd * two_electron_energy_ + (1.0 - lmbd) * (coulomb_energy + hartree_ex_energy + lrc_ex_energy) + mcpdft_xc_energy;
    // double wf_contribution  = one_electron_energy + lr_Vee_energy_ + mcpdft_lambda * sr_Vee_energy_;
    double wf_contribution  = one_electron_energy + mcpdft_lambda * two_electron_energy_;
    double dft_contribution = (1.0 - mcpdft_lambda) * (coulomb_energy_ + mcpdft_ex) + (1.0 - mcpdft_lambda * mcpdft_lambda) * mcpdft_ec;
    double total_energy = wf_contribution + dft_contribution + nuclear_repulsion_energy;
    if( options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT" || options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT" ) 
      total_energy += (mcpdft_lambda * hf_ex_energy_) + (mcpdft_lambda * mcpdft_lambda * mp2_corr_energy_);

    if( options_.get_str("MCPDFT_METHOD") == "MCPDFT") {
      outfile->Printf("    * MCPDFT total energy   =     ");
    }else if( options_.get_str("MCPDFT_METHOD") == "1H_MCPDFT") {
            outfile->Printf("    * 1H-MCPDFT total energy      =      ");
    }else if( options_.get_str("MCPDFT_METHOD") == "1DH_MCPDFT") {
            outfile->Printf("    * 1DH-MCPDFT total energy     =      ");
    }else if( options_.get_str("MCPDFT_METHOD") == "RS_MCPDFT") {
            outfile->Printf("    * RS-MCPDFT total energy      =      ");
    }else if( options_.get_str("MCPDFT_METHOD") == "RS1H_MCPDFT") {
            outfile->Printf("    * RS1H-MCPDFT total energy    =      ");
    }else if( options_.get_str("MCPDFT_METHOD") == "RS1DH_MCPDFT") {
            outfile->Printf("    * RS1DH-MCPDFT total energy   =      ");
    } 
    outfile->Printf("   %20.12lf\n\n",total_energy);

    return total_energy;
}

void MCPDFTSolver::BuildPi(double * D2ab) {

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    double ** phi = super_phi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (int mu = 0; mu < nmo_; mu++) {
            for (int nu = 0; nu < nmo_; nu++) {
                for (int lambda = 0; lambda < nmo_; lambda++) {
                    for (int sigma = 0; sigma < nmo_; sigma++) {

                        dum += phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];
                    }
                }
            }
        }

        pi_p[p] = dum;
    }

    if ( is_gga_ || is_meta_ ) {

       pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       
       double * pi_xp = pi_x_->pointer();
       double * pi_yp = pi_y_->pointer();
       double * pi_zp = pi_z_->pointer();
       
       double ** phi_x = super_phi_x_->pointer();
       double ** phi_y = super_phi_y_->pointer();
       double ** phi_z = super_phi_z_->pointer();

       for (int p = 0; p < phi_points_; p++) {

           double dum_x = 0.0;
           double dum_y = 0.0;
           double dum_z = 0.0;

           // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
           for (int mu = 0; mu < nmo_; mu++) {
               for (int nu = 0; nu < nmo_; nu++) {
                   for (int lambda = 0; lambda < nmo_; lambda++) {
                       for (int sigma = 0; sigma < nmo_; sigma++) {

                           dum_x += ( phi_x[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_x[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_x[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_x[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];

                           dum_y += ( phi_y[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_y[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_y[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_y[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];


                           dum_z += ( phi_z[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi_z[p][lambda] * phi[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi_z[p][sigma] * phi[p][nu] +
                                      phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi_z[p][nu] ) * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];
                       }
                   }
               }
           }

           pi_xp[p] = dum_x;
           pi_yp[p] = dum_y;
           pi_zp[p] = dum_z;

           // outfile->Printf("pi_x %15.15lf\n",pi_xp[p]);
           // outfile->Printf("pi_y %15.15lf\n",pi_yp[p]);
           // outfile->Printf("pi_z %15.15lf\n",pi_zp[p]);
           // outfile->Printf("\n    p");
           // outfile->Printf("    x[p]");
           // outfile->Printf("    y[p]");
           // outfile->Printf("    z[p]");
           // outfile->Printf("    pi\n\n");

           // for (int p = 0; p < phi_points_; p++) {

           //     outfile->Printf("    %d %20.5lf %20.5lf %20.5lf %20.15lf\n",p, grid_x_->pointer()[p], grid_y_->pointer()[p], grid_z_->pointer()[p], pi_p[p]);
           // }
       }
    }
}
void MCPDFTSolver::BuildPiLowMemory(tpdm * D2ab, int nab) {

    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    int nblocks = potential_->nblocks();

    Dimension phi_points_list = Dimension(nirrep_,"phi_points");
    for (int h = 0; h < nirrep_; h++) {
        phi_points_list[h] = max_points_;
    }
    std::shared_ptr<Matrix> myPhiAO(new Matrix(max_points_,nso_));
    std::shared_ptr<Matrix> myPhiMO(new Matrix(phi_points_list,nsopi_));
    std::shared_ptr<Matrix> tempPhi(new Matrix(phi_points_list,nsopi_));

    std::shared_ptr<Matrix> myPhiMO_x;
    std::shared_ptr<Matrix> myPhiMO_y;
    std::shared_ptr<Matrix> myPhiMO_z;

    // allocate memory for Pi, Pi_x, etc.
    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    double * pi_p  = pi_->pointer();

    double * pi_xp;
    double * pi_yp;
    double * pi_zp;

    if ( is_gga_ || is_meta_ ) {

        myPhiMO_x = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));
        myPhiMO_y = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));
        myPhiMO_z = (std::shared_ptr<Matrix>)(new Matrix(phi_points_list,nsopi_));

        pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        
        pi_xp = pi_x_->pointer();
        pi_yp = pi_y_->pointer();
        pi_zp = pi_z_->pointer();
    }

    // memory for rho
    rho_a_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p  = rho_->pointer();

    double * sigma_aap;
    double * sigma_bbp;
    double * sigma_abp;

    double * rho_a_xp;
    double * rho_b_xp;

    double * rho_a_yp;
    double * rho_b_yp;

    double * rho_a_zp;
    double * rho_b_zp;

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        sigma_aap = sigma_aa_->pointer();
        sigma_bbp = sigma_bb_->pointer();
        sigma_abp = sigma_ab_->pointer();

        rho_a_xp = rho_a_x_->pointer();
        rho_b_xp = rho_b_x_->pointer();

        rho_a_yp = rho_a_y_->pointer();
        rho_b_yp = rho_b_y_->pointer();

        rho_a_zp = rho_a_z_->pointer();
        rho_b_zp = rho_b_z_->pointer();

    }

    // number of nonzero elements in Da and Db:
    long int nb;
    long int na;

    std::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));
    psio->close(PSIF_V2RDM_D1A,1);

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));
    psio->close(PSIF_V2RDM_D1B,1);
        
    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential_->get_block(myblock);
        points_func_->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi   = points_func_->basis_value("PHI")->pointer();

        for (int p = 0; p < npoints; p++) {

            // grab block of AO-basis phi matrix
            myPhiAO->zero();
            for (int nu = 0; nu < nlocal; nu++) {
                int nug = function_map[nu];
                myPhiAO->pointer()[p][nug] = phi[p][nu];
            }

            // transform AO-basis block of Phi matrix to MO basis
            for (int h = 0; h < nirrep_; h++) {
                if (nsopi_[h] != 0 ) {
                    F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                }
            }

            for (int h = 0; h < nirrep_; h++) {
                if (nsopi_[h] != 0 ) {
                    F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO->pointer(h)[0],nsopi_[h]);
                }
            }

            // rho_a(r)
            double duma = 0.0;
            for (int n = 0; n < na; n++) {
                
                int i = opdm_a_[n].i;
                int j = opdm_a_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                duma += myPhiMO->pointer(hi)[p][ii] * 
                        myPhiMO->pointer(hj)[p][jj] * opdm_a_[n].val;
            
            }
            rho_ap[phi_points_ + p] = duma;

            // rho_b(r)
            double dumb = 0.0;
            for (int n = 0; n < nb; n++) {

                int i = opdm_b_[n].i;
                int j = opdm_b_[n].j;

                int hi = symmetry_[i];
                int hj = symmetry_[j];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];

                dumb += myPhiMO->pointer(hi)[p][ii] *
                        myPhiMO->pointer(hj)[p][jj] * opdm_b_[n].val;

            }
            rho_bp[phi_points_ + p] = dumb;

            rho_p[phi_points_ + p] = rho_ap[phi_points_ + p] + rho_bp[phi_points_ + p];


            // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
            double dum = 0.0;
            for (int n = 0; n < nab; n++) {

                int i = D2ab[n].i;
                int j = D2ab[n].j;
                int k = D2ab[n].k;
                int l = D2ab[n].l;

                int hi = symmetry_[i];
                int hj = symmetry_[j];
                int hk = symmetry_[k];
                int hl = symmetry_[l];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                int kk = k - pitzer_offset_[hk];
                int ll = l - pitzer_offset_[hl];

                dum += myPhiMO->pointer(hi)[p][ii] *
                       myPhiMO->pointer(hj)[p][jj] *
                       myPhiMO->pointer(hk)[p][kk] *
                       myPhiMO->pointer(hl)[p][ll] * D2ab[n].val;

            }

            pi_p[phi_points_ + p] = dum;

            if ( is_gga_ || is_meta_ ) {

                double ** phi_x = points_func_->basis_value("PHI_X")->pointer();
                double ** phi_y = points_func_->basis_value("PHI_Y")->pointer();
                double ** phi_z = points_func_->basis_value("PHI_Z")->pointer();

                // grab block of AO-basis phi_x matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_x[p][nu];
                }

                // transform AO-basis block of Phi_x matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_x->pointer(h)[0],nsopi_[h]);
                    }
                }

                // grab block of AO-basis phi_y matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_y[p][nu];
                }

                // transform AO-basis block of Phi_y matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_y->pointer(h)[0],nsopi_[h]);
                    }
                }
                // grab block of AO-basis phi_z matrix
                myPhiAO->zero();
                for (int nu = 0; nu < nlocal; nu++) {
                    int nug = function_map[nu];
                    myPhiAO->pointer()[p][nug] = phi_z[p][nu];
                }

                // transform AO-basis block of Phi_z matrix to MO basis
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nso_,1.0,ao2so->pointer(h)[0],nsopi_[h],myPhiAO->pointer()[0],nso_,0.0,tempPhi->pointer(h)[0],nsopi_[h]);
                    }
                }

                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_[h] != 0 ) {
                        F_DGEMM('n','n',nsopi_[h],npoints,nsopi_[h],1.0,Ca_->pointer(h)[0],nsopi_[h],tempPhi->pointer(h)[0],nsopi_[h],0.0,myPhiMO_z->pointer(h)[0],nsopi_[h]);
                    }
                }

                // rho'_a(r)

                double duma_x = 0.0;
                double duma_y = 0.0;
                double duma_z = 0.0;

                for (int n = 0; n < na; n++) {

                    int i = opdm_a_[n].i;
                    int j = opdm_a_[n].j;

                    int hi = symmetry_[i];
                    int hj = symmetry_[j];

                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];

                    duma_x += ( myPhiMO_x->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_x->pointer(hj)[p][jj] ) * opdm_a_[n].val;

                    duma_y += ( myPhiMO_y->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_y->pointer(hj)[p][jj] ) * opdm_a_[n].val;

                    duma_z += ( myPhiMO_z->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_z->pointer(hj)[p][jj] ) * opdm_a_[n].val;

                }

                // rho'_b(r)

                double dumb_x = 0.0;
                double dumb_y = 0.0;
                double dumb_z = 0.0;

                for (int n = 0; n < nb; n++) {

                    int i = opdm_b_[n].i;
                    int j = opdm_b_[n].j;

                    int hi = symmetry_[i];
                    int hj = symmetry_[j];

                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];

                    dumb_x += ( myPhiMO_x->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_x->pointer(hj)[p][jj] ) * opdm_b_[n].val;

                    dumb_y += ( myPhiMO_y->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_y->pointer(hj)[p][jj] ) * opdm_b_[n].val;

                    dumb_z += ( myPhiMO_z->pointer(hi)[p][ii] * myPhiMO->pointer(hj)[p][jj]
                            +   myPhiMO->pointer(hi)[p][ii] * myPhiMO_z->pointer(hj)[p][jj] ) * opdm_b_[n].val;

                }

                rho_a_xp[phi_points_ + p] = duma_x;
                rho_b_xp[phi_points_ + p] = dumb_x;

                rho_a_yp[phi_points_ + p] = duma_y;
                rho_b_yp[phi_points_ + p] = dumb_y;

                rho_a_zp[phi_points_ + p] = duma_z;
                rho_b_zp[phi_points_ + p] = dumb_z;

                sigma_aap[phi_points_ + p] = ( rho_a_xp[phi_points_ + p] * rho_a_xp[phi_points_ + p] ) 
                                           + ( rho_a_yp[phi_points_ + p] * rho_a_yp[phi_points_ + p] )
                                           + ( rho_a_zp[phi_points_ + p] * rho_a_zp[phi_points_ + p] );
                sigma_bbp[phi_points_ + p] = ( rho_b_xp[phi_points_ + p] * rho_b_xp[phi_points_ + p] ) 
                                           + ( rho_b_yp[phi_points_ + p] * rho_b_yp[phi_points_ + p] ) 
                                           + ( rho_b_zp[phi_points_ + p] * rho_b_zp[phi_points_ + p] );
                sigma_abp[phi_points_ + p] = ( rho_a_xp[phi_points_ + p] * rho_b_xp[phi_points_ + p] ) 
                                           + ( rho_a_yp[phi_points_ + p] * rho_b_yp[phi_points_ + p] ) 
                                           + ( rho_a_zp[phi_points_ + p] * rho_b_zp[phi_points_ + p] );

                double dum_x = 0.0;
                double dum_y = 0.0;
                double dum_z = 0.0;

                // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
                for (int n = 0; n < nab; n++) {

                    int i = D2ab[n].i;
                    int j = D2ab[n].j;
                    int k = D2ab[n].k;
                    int l = D2ab[n].l;
                    
                    int hi = symmetry_[i];
                    int hj = symmetry_[j];
                    int hk = symmetry_[k]; 
                    int hl = symmetry_[l];
                    
                    int ii = i - pitzer_offset_[hi];
                    int jj = j - pitzer_offset_[hj];
                    int kk = k - pitzer_offset_[hk];
                    int ll = l - pitzer_offset_[hl];

                    dum_x += ( myPhiMO_x->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_x->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_x->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_x->pointer(hl)[p][ll] ) * D2ab[n].val;

                    dum_y += ( myPhiMO_y->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_y->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_y->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_y->pointer(hl)[p][ll] ) * D2ab[n].val;

                    dum_z += ( myPhiMO_z->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO_z->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO_z->pointer(hk)[p][kk] * 
                               myPhiMO->pointer(hl)[p][ll] +

                               myPhiMO->pointer(hi)[p][ii] *
                               myPhiMO->pointer(hj)[p][jj] *
                               myPhiMO->pointer(hk)[p][kk] * 
                               myPhiMO_z->pointer(hl)[p][ll] ) * D2ab[n].val;
                }

                pi_xp[phi_points_ + p] = dum_x;
                pi_yp[phi_points_ + p] = dum_y;
                pi_zp[phi_points_ + p] = dum_z;
            }

        }
        phi_points_ += npoints;
    }

    free(opdm_a_);
    free(opdm_b_);
}

void MCPDFTSolver::BuildExchangeCorrelationHole(int p, tpdm * D2ab, int nab, tpdm * D2aa, int naa, tpdm * D2bb, int nbb) {

    double integral = 0.0;
    double * rho_p = rho_->pointer();
    double * w_p   = grid_w_->pointer();
    for (int q = 0; q < phi_points_; q++) {

        // pi(r,r') = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r',nu) * phi(r,lambda) * phi(r',sigma)
        double pi = 0.0;
        for (int n = 0; n < nab; n++) {

            int i = D2ab[n].i;
            int j = D2ab[n].j;
            int k = D2ab[n].k;
            int l = D2ab[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * D2ab[n].val * 0.5;

            pi += super_phi_->pointer(hi)[q][ii] * 
                  super_phi_->pointer(hj)[p][jj] * 
                  super_phi_->pointer(hk)[q][kk] * 
                  super_phi_->pointer(hl)[p][ll] * D2ab[n].val * 0.5;

        }
        for (int n = 0; n < naa; n++) {

            int i = D2aa[n].i;
            int j = D2aa[n].j;
            int k = D2aa[n].k;
            int l = D2aa[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * D2aa[n].val * 0.5;

        }
        for (int n = 0; n < nbb; n++) {

            int i = D2bb[n].i;
            int j = D2bb[n].j;
            int k = D2bb[n].k;
            int l = D2bb[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            pi += super_phi_->pointer(hi)[p][ii] * 
                  super_phi_->pointer(hj)[q][jj] * 
                  super_phi_->pointer(hk)[p][kk] * 
                  super_phi_->pointer(hl)[q][ll] * D2bb[n].val * 0.5;

        }
        double nxc = (2.0 * pi - rho_p[p] * rho_p[q]) / rho_p[p];
        //if ( fabs(grid_x_->pointer()[q] < 1e-6 ) ) {
        //printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",grid_x_->pointer()[q],grid_y_->pointer()[q],grid_z_->pointer()[q],nxc);
        //}

        integral += nxc * w_p[q];

    }
    printf("%20.12lf\n",integral);
    exit(0);

}

void MCPDFTSolver::BuildPiFast(tpdm * D2ab, int nab) {

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (int n = 0; n < nab; n++) {

            int i = D2ab[n].i;
            int j = D2ab[n].j;
            int k = D2ab[n].k;
            int l = D2ab[n].l;

            int hi = symmetry_[i];
            int hj = symmetry_[j];
            int hk = symmetry_[k];
            int hl = symmetry_[l];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];
            int kk = k - pitzer_offset_[hk];
            int ll = l - pitzer_offset_[hl];

            dum += super_phi_->pointer(hi)[p][ii] * 
                   super_phi_->pointer(hj)[p][jj] * 
                   super_phi_->pointer(hk)[p][kk] * 
                   super_phi_->pointer(hl)[p][ll] * D2ab[n].val;

        }

        pi_p[p] = dum;
    }

    if ( is_gga_ || is_meta_ ) {

        pi_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        pi_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        
        double * pi_xp = pi_x_->pointer();
        double * pi_yp = pi_y_->pointer();
        double * pi_zp = pi_z_->pointer();
        
        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;

            // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
            for (int n = 0; n < nab; n++) {

                int i = D2ab[n].i;
                int j = D2ab[n].j;
                int k = D2ab[n].k;
                int l = D2ab[n].l;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                int hk = symmetry_[k]; 
                int hl = symmetry_[l];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                int kk = k - pitzer_offset_[hk];
                int ll = l - pitzer_offset_[hl];

                dum_x += ( super_phi_x_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_x_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_x_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_x_->pointer(hl)[p][ll] ) * D2ab[n].val;

                dum_y += ( super_phi_y_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_y_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_y_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_y_->pointer(hl)[p][ll] ) * D2ab[n].val;

                dum_z += ( super_phi_z_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_z_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_z_->pointer(hk)[p][kk] * 
                           super_phi_->pointer(hl)[p][ll] +

                           super_phi_->pointer(hi)[p][ii] *
                           super_phi_->pointer(hj)[p][jj] *
                           super_phi_->pointer(hk)[p][kk] * 
                           super_phi_z_->pointer(hl)[p][ll] ) * D2ab[n].val;

            }

            pi_xp[p] = dum_x;
            pi_yp[p] = dum_y;
            pi_zp[p] = dum_z;

            // outfile->Printf("pi_x %15.15lf\n",pi_xp[p]);
            // outfile->Printf("pi_y %15.15lf\n",pi_yp[p]);
            // outfile->Printf("pi_z %15.15lf\n",pi_zp[p]);
            // outfile->Printf("\n    p");
            // outfile->Printf("    x[p]");
            // outfile->Printf("    y[p]");
            // outfile->Printf("    z[p]");
            // outfile->Printf("    pi\n\n");

            // for (int p = 0; p < phi_points_; p++) {

            //     outfile->Printf("    %d %20.5lf %20.5lf %20.5lf %20.15lf\n",p, grid_x_->pointer()[p], grid_y_->pointer()[p], grid_z_->pointer()[p], pi_p[p]);
            // }
        }
    }
}

void MCPDFTSolver::BuildRho() {

    rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // m_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    // zeta_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // rs_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double ** phi   = super_phi_->pointer();

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p = rho_->pointer();
    // double * m_p = m_->pointer();
    // 
    // double * zeta_p = zeta_->pointer();
    // double * rs_p = rs_->pointer();
    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;

    double ** dap = Da_->pointer();
    double ** dbp = Db_->pointer();

    for (int p = 0; p < phi_points_; p++) {
        double duma   = 0.0;
        double dumb   = 0.0;
        for (int sigma = 0; sigma < nmo_; sigma++) {
            for (int nu = 0; nu < nmo_; nu++) {
                duma += phi[p][sigma] * phi[p][nu] * dap[sigma][nu];
                dumb += phi[p][sigma] * phi[p][nu] * dbp[sigma][nu];
            }
        }
        rho_ap[p] = duma;
        rho_bp[p] = dumb;
        rho_p[p] = rho_ap[p] + rho_bp[p];    

        temp_tot += rho_p[p] * grid_w_->pointer()[p];
        temp_a += rho_ap[p] * grid_w_->pointer()[p];
        temp_b += rho_bp[p] * grid_w_->pointer()[p];

    }
    outfile->Printf("\n");
    outfile->Printf("      Integrated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        // tau_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        // tau_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        double * sigma_aap = sigma_aa_->pointer();
        double * sigma_bbp = sigma_bb_->pointer();
        double * sigma_abp = sigma_ab_->pointer();
        
        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();

        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        // double * tau_ap = tau_a_->pointer();
        // double * tau_bp = tau_b_->pointer();

        for (int p = 0; p < phi_points_; p++) {
            double duma_x = 0.0;
            double dumb_x = 0.0;
            double duma_y = 0.0;
            double dumb_y = 0.0;
            double duma_z = 0.0;
            double dumb_z = 0.0;
            // double dumta = 0.0;
            // double dumtb = 0.0;

            for (int sigma = 0; sigma < nmo_; sigma++) {
                for (int nu = 0; nu < nmo_; nu++) {
                    duma_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * dap[sigma][nu];
                    dumb_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * dbp[sigma][nu];

                    duma_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * dap[sigma][nu];
                    dumb_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * dbp[sigma][nu];

                    duma_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * dap[sigma][nu];
                    dumb_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * dbp[sigma][nu];
                }
            }
            rho_a_xp[p] = duma_x;
            rho_b_xp[p] = dumb_x;

            rho_a_yp[p] = duma_y;
            rho_b_yp[p] = dumb_y;

            rho_a_zp[p] = duma_z;
            rho_b_zp[p] = dumb_z;

            sigma_aap[p] = ( rho_a_xp[p] * rho_a_xp[p] ) +  ( rho_a_yp[p] * rho_a_yp[p] ) + ( rho_a_zp[p] * rho_a_zp[p] );
            sigma_bbp[p] = ( rho_b_xp[p] * rho_b_xp[p] ) +  ( rho_b_yp[p] * rho_b_yp[p] ) + ( rho_b_zp[p] * rho_b_zp[p] );
            sigma_abp[p] = ( rho_a_xp[p] * rho_b_xp[p] ) +  ( rho_a_yp[p] * rho_b_yp[p] ) + ( rho_a_zp[p] * rho_b_zp[p] );

            // tau_ap[p] = dumta;
            // tau_bp[p] = dumtb;

        }
    }
}

void MCPDFTSolver::BuildRhoFast(int na, int nb) {

    rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // m_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    // zeta_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // rs_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double ** phi   = super_phi_->pointer();

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * rho_p = rho_->pointer();
    // double * m_p = m_->pointer();
    // 
    // double * zeta_p = zeta_->pointer();
    // double * rs_p = rs_->pointer();
    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        // rho_a(r)
        double duma = 0.0;
        for (int n = 0; n < na; n++) {

            int i = opdm_a_[n].i;
            int j = opdm_a_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            duma += super_phi_->pointer(hi)[p][ii] *
                    super_phi_->pointer(hj)[p][jj] * opdm_a_[n].val;

        }
        rho_ap[p] = duma;

        // rho_b(r)
        double dumb = 0.0;
        for (int n = 0; n < nb; n++) {

            int i = opdm_b_[n].i;
            int j = opdm_b_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            dumb += super_phi_->pointer(hi)[p][ii] *
                    super_phi_->pointer(hj)[p][jj] * opdm_b_[n].val;

        }
        rho_bp[p] = dumb;

        rho_p[p] = rho_ap[p] + rho_bp[p];    

        temp_tot += rho_p[p]  * grid_w_->pointer()[p];
        temp_a   += rho_ap[p] * grid_w_->pointer()[p];
        temp_b   += rho_bp[p] * grid_w_->pointer()[p];

        // rho_p[p] = rho_ap[p] + rho_bp[p];
        // m_p[p] =  rho_ap[p] - rho_bp[p];

        // zeta_p[p] =  m_p[p]  / rho_p[p];
        // rs_p[p] = pow( 3.0 / ( 4.0 * M_PI * rho_p[p] ) , 1.0/3.0 );
    }
    outfile->Printf("\n");
    outfile->Printf("      Integrated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        // tau_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        // tau_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        double * sigma_aap = sigma_aa_->pointer();
        double * sigma_bbp = sigma_bb_->pointer();
        double * sigma_abp = sigma_ab_->pointer();

        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();

        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        // double * tau_ap = tau_a_->pointer();
        // double * tau_bp = tau_b_->pointer();

        for (int p = 0; p < phi_points_; p++) {

            // rho'_a(r)

            double duma_x = 0.0;
            double duma_y = 0.0;
            double duma_z = 0.0;

            for (int n = 0; n < na; n++) {

                int i = opdm_a_[n].i;
                int j = opdm_a_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                duma_x += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] ) * opdm_a_[n].val;

                duma_y += ( super_phi_y_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj] ) * opdm_a_[n].val;

                duma_z += ( super_phi_z_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_a_[n].val;

            }

            // rho'_b(r)

            double dumb_x = 0.0;
            double dumb_y = 0.0;
            double dumb_z = 0.0;

            for (int n = 0; n < nb; n++) {
                
                int i = opdm_b_[n].i;
                int j = opdm_b_[n].j;
                
                int hi = symmetry_[i];
                int hj = symmetry_[j];
                
                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];
                
                dumb_x += ( super_phi_x_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_x_->pointer(hj)[p][jj] ) * opdm_b_[n].val;
                
                dumb_y += ( super_phi_y_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_y_->pointer(hj)[p][jj] ) * opdm_b_[n].val;
                
                dumb_z += ( super_phi_z_->pointer(hi)[p][ii] * super_phi_->pointer(hj)[p][jj] 
                        +   super_phi_->pointer(hi)[p][ii] * super_phi_z_->pointer(hj)[p][jj] ) * opdm_b_[n].val;

            }

            rho_a_xp[p] = duma_x;
            rho_b_xp[p] = dumb_x;

            rho_a_yp[p] = duma_y;
            rho_b_yp[p] = dumb_y;

            rho_a_zp[p] = duma_z;
            rho_b_zp[p] = dumb_z;

            sigma_aap[p] = ( rho_a_xp[p] * rho_a_xp[p] ) +  ( rho_a_yp[p] * rho_a_yp[p] ) + ( rho_a_zp[p] * rho_a_zp[p] );
            sigma_bbp[p] = ( rho_b_xp[p] * rho_b_xp[p] ) +  ( rho_b_yp[p] * rho_b_yp[p] ) + ( rho_b_zp[p] * rho_b_zp[p] );
            sigma_abp[p] = ( rho_a_xp[p] * rho_b_xp[p] ) +  ( rho_a_yp[p] * rho_b_yp[p] ) + ( rho_a_zp[p] * rho_b_zp[p] );

            // tau_ap[p] = dumta;
            // tau_bp[p] = dumtb;

        }
    }
}

double MCPDFTSolver::RangeSeparatedTEE(std::string range_separation_type) {

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);

    //Ca_->print();

    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    // ab
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    long int nab;
    psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

    for (int n = 0; n < nab; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2ab[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AB,1);

    // aa
    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

    long int naa;
    psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

    for (int n = 0; n < naa; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2aa[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AA,1);

    // bb
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

    long int nbb;
    psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

    for (int n = 0; n < nbb; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2bb[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2BB,1);

/*
    SharedMatrix eri (new Matrix(mints->mo_erfc_eri(options_.get_double("MCPDFT_OMEGA"),Ca_,Cb_,Ca_,Cb_)));
    double ** eri_p = eri->pointer();

    // (11|22)
    // (ij|kl) = erf_eri->pointer()[i*nmo_+j][k*nmo_+l]
    double e2 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            int ij = i*nmo_+j;
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                    int kl = k*nmo_+l;
                    e2 += 0.5 * eri_p[ij][kl] * D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                    e2 += 0.5 * eri_p[ij][kl] * D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                    e2 += eri_p[ij][kl] * D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                }
            }
        }
    }
*/
    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    if (range_separation_type == "LR") {
       outfile->Printf("    ==> Transform ERF integrals <==\n");
       outfile->Printf("\n");
       // write erf integrals to disk
       mints->integrals_erf(options_.get_double("MCPDFT_OMEGA"));
    }else if (range_separation_type == "SR") {
             outfile->Printf("    ==> Transform ERFC integrals <==\n");
             outfile->Printf("\n");
             // write erfc integrals to disk
             mints->integrals_erfc(options_.get_double("MCPDFT_OMEGA"));
    }else{
         throw PsiException("The argument of RangeSeparatedTEI() function should either be \"SR\" or \"LR\"",__FILE__,__LINE__);
    } 

    
    // transform range_separated (erf/erfc) integrals
    double start = omp_get_wtime();

    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);

    std::shared_ptr<IntegralTransform> ints(new IntegralTransform(reference_wavefunction_, spaces,
        IntegralTransform::TransformationType::Restricted, IntegralTransform::OutputType::IWLOnly, 
        IntegralTransform::MOOrdering::PitzerOrder, IntegralTransform::FrozenOrbitals::None, false));

    ints->set_dpd_id(0);
    ints->set_keep_iwl_so_ints(true);
    ints->set_keep_dpd_so_ints(true);
    if (range_separation_type == "LR") {
       ints->set_so_tei_file(PSIF_SO_ERF_TEI);
    }else if (range_separation_type == "SR") {
             ints->set_so_tei_file(PSIF_SO_ERFC_TEI);
    }else{
         throw PsiException("The argument of RangeSeparatedTEI() function should either be \"SR\" or \"LR\"",__FILE__,__LINE__);
    } 

    ints->initialize();

    ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

    double end = omp_get_wtime();

    outfile->Printf("\n");
    outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
    outfile->Printf("\n");

    // read range-separated integrals from disk
    ReadRangeSeparatedIntegrals();

    double range_separated_2e_energy = 0.0;
    for (int i = 0; i < nmo_; i++) {

        int hi = symmetry_[i];

        for (int j = 0; j < nmo_; j++) {

            int hj  = symmetry_[j];
            int hij = hi ^ hj;
            int ij  = ibas_[hij][i][j];

            for (int k = 0; k < nmo_; k++) {

                int hk = symmetry_[k];

                for (int l = 0; l < nmo_; l++) {

                    int hl  = symmetry_[l];
                    int hkl = hk ^ hl;

                    if ( hij != hkl ) continue;

                    int kl = ibas_[hkl][k][l];

                    long int offset = 0;
                    for (int myh = 0; myh < hij; myh++) {
                        offset += (long int)gems_[myh].size() * ( (long int)gems_[myh].size() + 1 ) / 2;
                    }
                    double val = erfc_tei_[offset + INDEX(ij,kl)];

                    range_separated_2e_energy += 0.5 * val * D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                    range_separated_2e_energy += 0.5 * val * D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                    range_separated_2e_energy +=       val * D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+l];
                }
            }
        }
    }
    // printf("  erfc coulomb energy: %20.12lf %20.12lf\n",erfc_coulomb_energy, two_electron_energy_);

    free(D2aa);
    free(D2ab);
    free(D2bb);
 
    return range_separated_2e_energy;
}

std::vector< std::shared_ptr<Matrix> > MCPDFTSolver::BuildJK() {

    std::shared_ptr<PSIO> psio (new PSIO());

    // psio->set_pid("18332");

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm_a_ = (opdm *)malloc(na * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a_,na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    for (int n = 0; n < na; n++) {

        int i = opdm_a_[n].i;
        int j = opdm_a_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the alpha OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Da_->pointer(hi)[ii][jj] = opdm_a_[n].val;

    }

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm_b_ = (opdm *)malloc(nb * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b_,nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

   for (int n = 0; n < nb; n++) {

        int i = opdm_b_[n].i;
        int j = opdm_b_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Db_->pointer(hi)[ii][jj] = opdm_b_[n].val;

    }

    free(opdm_a_);
    free(opdm_b_);

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    std::vector< std::shared_ptr<Matrix> > JK;

    // JK object (note this is hard-coded to use density fitting ...)
    if (options_.get_str("MCPDFT_TYPE") == "DF") {

         // get auxiliary basis:
         std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");
        
        // outfile->Printf("\n");
        // outfile->Printf("    ==> The JK type for the MCPDFT calculation is: %s", options_.get_str("MCPDFT_TYPE"));
        // outfile->Printf(" <==\n");

        std::shared_ptr<DiskDFJK> jk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary));
        
        // memory for jk (say, 85% of what is available)
        jk->set_memory(0.90 * Process::environment.get_memory());
        // jk->set_memory(0.85 * available_memory_);

        // integral cutoff
        jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        jk->set_do_J(true);
        jk->set_do_K(false);
        jk->set_do_wK(false);
        jk->set_omega(false);

        if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {

           if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
              && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {

              jk->set_do_K(true);
              jk->set_do_wK(true);
              jk->set_omega(options_.get_double("MCPDFT_OMEGA"));
              // scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
              // std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
              // jk->set_omega(functional->x_omega());
              jk->set_do_K(true);
           }
           jk->set_do_K(true);
        }

        jk->initialize();

        std::vector<SharedMatrix>& C_left  = jk->C_left();
        std::vector<SharedMatrix>& C_right = jk->C_right();

        // use alpha and beta C matrices

        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );

        C_left.clear();
        C_right.clear();

// std::shared_ptr<Matrix> D (new Matrix(Da_));
// D->add(Db_);


        // alpha first 
        myCa->zero();
        myCa->gemm('t','n',1.0,Da_,Ca_,0.0);
        myCa->transpose_this();
        C_left.push_back(myCa);
        C_right.push_back(Ca_);

        // beta second 
        myCb->zero();
        myCb->gemm('t','n',1.0,Db_,Cb_,0.0);
        myCb->transpose_this();
        C_left.push_back(myCb);
        C_right.push_back(Cb_);

        // Let jk compute for the given C_left/C_right

        jk->compute();

        std::shared_ptr<Matrix> Ja = jk->J()[0];
        std::shared_ptr<Matrix> Jb = jk->J()[1];

        Ja->transform(Ca_);
        Jb->transform(Cb_);

        JK.push_back(Ja);
        JK.push_back(Jb);

        if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {

            std::shared_ptr<Matrix> Ka = jk->K()[0];
            std::shared_ptr<Matrix> Kb = jk->K()[1];
            
            Ka->transform(Ca_);
            Kb->transform(Cb_);
            
            JK.push_back(Ka);
            JK.push_back(Kb);

            if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
               && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {
    
               std::shared_ptr<Matrix> wKa = jk->wK()[0];
               std::shared_ptr<Matrix> wKb = jk->wK()[1];
               
               wKa->transform(Ca_);
               wKb->transform(Cb_);
               
               JK.push_back(wKa);
               JK.push_back(wKb);
            }
        }
        jk.reset();
        return JK;

    }else if (options_.get_str("MCPDFT_TYPE") == "PK") {

        // outfile->Printf("\n");
        // outfile->Printf("    ==> The JK type for the MCPDFT calculation is: %s", options_.get_str("MCPDFT_TYPE"));
        // outfile->Printf(" <==\n");

        std::shared_ptr<PKJK> jk = (std::shared_ptr<PKJK>)(new PKJK(primary,options_));

        // memory for jk (say, 85% of what is available)
        jk->set_memory(0.90 * Process::environment.get_memory());
        // jk->set_memory(0.85 * available_memory_);

        // integral cutoff
        jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        jk->set_do_J(true);
        jk->set_do_K(false);
        jk->set_do_wK(false);
        jk->set_omega(false);

        if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {

           if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
              && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {

              scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
              std::shared_ptr<SuperFunctional> functional = scfwfn->functional();

              jk->set_do_K(true);
              jk->set_do_wK(true);
              jk->set_omega(options_.get_double("MCPDFT_OMEGA"));
              // jk->set_omega(functional->x_omega());
              jk->set_do_K(true);
           }
           jk->set_do_K(true);
        }

        jk->initialize();

        std::vector<SharedMatrix>& C_left  = jk->C_left();
        std::vector<SharedMatrix>& C_right = jk->C_right();

        // use alpha and beta C matrices

        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );

        C_left.clear();
        C_right.clear();

        // alpha first 
        myCa->zero();
        myCa->gemm('t','n',1.0,Da_,Ca_,0.0);
        myCa->transpose_this();
        C_left.push_back(myCa);
        C_right.push_back(Ca_);

        // beta second 
        myCb->zero();
        myCb->gemm('t','n',1.0,Db_,Cb_,0.0);
        myCb->transpose_this();
        C_left.push_back(myCb);
        C_right.push_back(Cb_);

        // Let jk compute for the given C_left/C_right

        jk->compute();

        std::shared_ptr<Matrix> Ja = jk->J()[0];
        std::shared_ptr<Matrix> Jb = jk->J()[1];

        Ja->transform(Ca_);
        Jb->transform(Cb_);

        JK.push_back(Ja);
        JK.push_back(Jb);

        if ( (options_.get_str("MCPDFT_METHOD") != "MCPDFT") ) {

           std::shared_ptr<Matrix> Ka = jk->K()[0];
           std::shared_ptr<Matrix> Kb = jk->K()[1];
           
           Ka->transform(Ca_);
           Kb->transform(Cb_);
           
           JK.push_back(Ka);
           JK.push_back(Kb);
        }

        if (  (options_.get_str("MCPDFT_METHOD") != "1H_MCPDFT") 
           && (options_.get_str("MCPDFT_METHOD") != "1DH_MCPDFT") ) {

           std::shared_ptr<Matrix> wKa = jk->wK()[0];
           std::shared_ptr<Matrix> wKb = jk->wK()[1];
           
           wKa->transform(Ca_);
           wKb->transform(Cb_);
           
           JK.push_back(wKa);
           JK.push_back(wKb);
        }

        return JK;
        
    }else {
        throw PsiException("invalid MCSCF_TYPE",__FILE__,__LINE__);
    }
}

void MCPDFTSolver::Build_R(){

    double tol = 1.0e-20;

    R_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    double * R_p = R_->pointer();
    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
 
    for (int p = 0; p < phi_points_; p++) {
        
        R_p[p] = 4 * pi_p[p] / (rho_p[p] * rho_p[p]); 
    }
} 

void MCPDFTSolver::Translate(){

    double tol = 1.0e-20;

    tr_rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * tr_rho_ap = tr_rho_a_->pointer();
    double * tr_rho_bp = tr_rho_b_->pointer();

    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();

    // if ( deriv_ != 0) {

    // tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    // tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    // tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    // 
    // double * rho_a_xp = rho_a_x_->pointer();
    // double * rho_b_xp = rho_b_x_->pointer();

    // double * rho_a_yp = rho_a_y_->pointer();
    // double * rho_b_yp = rho_b_y_->pointer();

    // double * rho_a_zp = rho_a_z_->pointer();
    // double * rho_b_zp = rho_b_z_->pointer();

    // double * pi_xp = pi_x_->pointer();
    // double * pi_yp = pi_y_->pointer();
    // double * pi_zp = pi_z_->pointer();

    // double * tr_rho_a_xp = tr_rho_a_x_->pointer();
    // double * tr_rho_b_xp = tr_rho_b_x_->pointer();

    // double * tr_rho_a_yp = tr_rho_a_y_->pointer();
    // double * tr_rho_b_yp = tr_rho_b_y_->pointer();

    // double * tr_rho_a_zp = tr_rho_a_z_->pointer();
    // double * tr_rho_b_zp = tr_rho_b_z_->pointer();

    // }

    double temp_tot = 0.0;
    double temp_spin = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_p[p];
        double pi = pi_p[p];
        
        double zeta = 0.0;
        double R = 0.0;

	// if ( !(rho < tol) && !(pi < tol) ) {
	if ( !(rho < tol) && !(pi < 0.0) ) {

           R = R_p[p];
           // R = tanh(R);
           // outfile->Printf("tanh(R) = %12.5lf\n",R);

           if ( (1.0 - R) > tol ) {

              zeta = sqrt(1.0 - R);

           }else{

                zeta = 0.0;
           }

           tr_rho_ap[p] = (1.0 + zeta) * (rho/2.0);
           tr_rho_bp[p] = (1.0 - zeta) * (rho/2.0);
           // tr_rs_p[p] = pow( 3.0 / ( 4.0 * M_PI * (tr_rho_ap[p] + tr_rho_bp[p]) ) , 1.0/3.0 );

        }else {

               tr_rho_ap[p] = 0.0;
               tr_rho_bp[p] = 0.0;
        }

        temp_tot  += (tr_rho_ap[p] + tr_rho_bp[p]) * grid_w_->pointer()[p];
        temp_spin += (tr_rho_ap[p] - tr_rho_bp[p]) * grid_w_->pointer()[p];
        temp_a += tr_rho_ap[p] * grid_w_->pointer()[p];
        temp_b += tr_rho_bp[p] * grid_w_->pointer()[p];

        // outfile->Printf("zeta = %12.15lf\n",tr_zeta_p[p]);
        // outfile->Printf("rs = %12.15lf\n",tr_rs_p[p]);
        
        // tr_m_p[p] = (R_p[p] > 1.0) ? 0.0 : rho_p[p] * sqrt(1.0 - R_p[p]);

        // tr_rho_ap[p] = ( (R_p[p] > 1.0) ? (rho/2.0) : (rho/2.0) * ( 1.0 + sqrt(1.0 - R_p[p]) ) );
        // tr_rho_bp[p] = ( (R_p[p] > 1.0) ? (rho/2.0) : (rho/2.0) * ( 1.0 - sqrt(1.0 - R_p[p]) ) );

        // tr_zeta_p[p] =  ( tr_rho_ap[p] - tr_rho_bp[p] ) / ( tr_rho_ap[p] + tr_rho_bp[p] );
        // tr_rs_p[p] = pow( 3.0 / ( 4.0 * M_PI * (tr_rho_ap[p] + tr_rho_bp[p]) ) , 1.0/3.0 );
    }

    outfile->Printf("\n");
    outfile->Printf("      Integrated translated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated translated spin density  = %20.12lf\n",temp_spin);
    outfile->Printf("      Integrated translated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated translated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

       tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       double * rho_a_xp = rho_a_x_->pointer();
       double * rho_b_xp = rho_b_x_->pointer();

       double * rho_a_yp = rho_a_y_->pointer();
       double * rho_b_yp = rho_b_y_->pointer();

       double * rho_a_zp = rho_a_z_->pointer();
       double * rho_b_zp = rho_b_z_->pointer();

       // double * pi_xp = pi_x_->pointer();
       // double * pi_yp = pi_y_->pointer();
       // double * pi_zp = pi_z_->pointer();

       double * tr_rho_a_xp = tr_rho_a_x_->pointer();
       double * tr_rho_b_xp = tr_rho_b_x_->pointer();

       double * tr_rho_a_yp = tr_rho_a_y_->pointer();
       double * tr_rho_b_yp = tr_rho_b_y_->pointer();

       double * tr_rho_a_zp = tr_rho_a_z_->pointer();
       double * tr_rho_b_zp = tr_rho_b_z_->pointer();

       double * tr_sigma_aap = tr_sigma_aa_->pointer();
       double * tr_sigma_abp = tr_sigma_ab_->pointer();
       double * tr_sigma_bbp = tr_sigma_bb_->pointer();

       for (int p = 0; p < phi_points_; p++) {

           double rho = rho_p[p];
           double pi = pi_p[p];

           double rho_x = rho_a_xp[p] + rho_b_xp[p];
           double rho_y = rho_a_yp[p] + rho_b_yp[p];
           double rho_z = rho_a_zp[p] + rho_b_zp[p];
           
           double zeta = 0.0;
           double R = 0.0;

           // if ( !(rho < tol) && !(pi < tol) ) {
	   if ( !(rho < tol) && !(pi < 0.0) ) {

              R = R_p[p];
              // R = tanh(R);
               
              if ( (1.0 - R) > tol )  {

                 zeta = sqrt(1.0 - R);

                 // tr_rho_a_xp[p] =  zeta_x * (rho/2.0) + (1.0 + zeta) * (rho_x/2.0);
                 // tr_rho_b_xp[p] = -zeta_x * (rho/2.0) + (1.0 - zeta) * (rho_x/2.0);

                 // tr_rho_a_yp[p] =  zeta_y * (rho/2.0) + (1.0 + zeta) * (rho_y/2.0);
                 // tr_rho_b_yp[p] = -zeta_y * (rho/2.0) + (1.0 - zeta) * (rho_y/2.0);

                 // tr_rho_a_zp[p] =  zeta_z * (rho/2.0) + (1.0 + zeta) * (rho_z/2.0);
                 // tr_rho_b_zp[p] = -zeta_z * (rho/2.0) + (1.0 - zeta) * (rho_z/2.0);

              }else{

                   zeta = 0.0;
              }

              tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0);
              tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0);

              tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0);
              tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0);

              tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0);
              tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0);

           }else {

                 tr_rho_a_xp[p] = 0.0;
                 tr_rho_b_xp[p] = 0.0;

                 tr_rho_a_yp[p] = 0.0;
                 tr_rho_b_yp[p] = 0.0;

                 tr_rho_a_zp[p] = 0.0;
                 tr_rho_b_zp[p] = 0.0;
           }

           tr_sigma_aap[p] = (tr_rho_a_xp[p] * tr_rho_a_xp[p]) + (tr_rho_a_yp[p] * tr_rho_a_yp[p]) + (tr_rho_a_zp[p] * tr_rho_a_zp[p]);
           tr_sigma_abp[p] = (tr_rho_a_xp[p] * tr_rho_b_xp[p]) + (tr_rho_a_yp[p] * tr_rho_b_yp[p]) + (tr_rho_a_zp[p] * tr_rho_b_zp[p]);
           tr_sigma_bbp[p] = (tr_rho_b_xp[p] * tr_rho_b_xp[p]) + (tr_rho_b_yp[p] * tr_rho_b_yp[p]) + (tr_rho_b_zp[p] * tr_rho_b_zp[p]);
       }
    }
    
}

void MCPDFTSolver::Fully_Translate(){
    
    double tol = 1.0e-20;

    double const R0 = 0.9;
    double const R1 = 1.15;
    double const A = -475.60656009;
    double const B = -379.47331922;
    double const C = -85.38149682;

    tr_rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    double * tr_rho_ap = tr_rho_a_->pointer();
    double * tr_rho_bp = tr_rho_b_->pointer();
    
    double * rho_p = rho_->pointer();
    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();

    double temp_tot = 0.0;
    double temp_a = 0.0;
    double temp_b = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_p[p];
        double pi = pi_p[p];
        double DelR = R_p[p] - R1;

        double zeta = 0.0;
        double R = 0.0;

	// if ( !(rho < tol) && !(pi < tol) ) {
	if ( !(rho < tol) && !(pi < 0.0) ) {

           R = R_p[p];
           // R = tanh(R);

           if ( ((1.0 - R) > tol) && ( R < R0 ) ) {

              zeta = sqrt(1.0 - R);

           }else if( !(R < R0) && !(R > R1) ) {

                   zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

           }else if( R > R1 ) {

                   zeta = 0.0;
           }

           tr_rho_ap[p] = (1.0 + zeta) * (rho/2.0);
           tr_rho_bp[p] = (1.0 - zeta) * (rho/2.0);

        }else{

             tr_rho_ap[p] = 0.0;
             tr_rho_bp[p] = 0.0;
        }

        temp_a += tr_rho_ap[p] * grid_w_->pointer()[p];
        temp_b += tr_rho_bp[p] * grid_w_->pointer()[p];
        temp_tot += ( tr_rho_bp[p] + tr_rho_ap[p] ) * grid_w_->pointer()[p];

    }

    outfile->Printf("\n");
    outfile->Printf("      Integrated fully translated total density = %20.12lf\n",temp_tot);
    outfile->Printf("      Integrated fully translated alpha density = %20.12lf\n",temp_a);
    outfile->Printf("      Integrated fully translated beta density  = %20.12lf\n",temp_b);
    outfile->Printf("\n");

    if ( is_gga_ || is_meta_ ) {

       tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       
       tr_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       double * tr_rho_a_xp = tr_rho_a_x_->pointer();
       double * tr_rho_b_xp = tr_rho_b_x_->pointer();

       double * tr_rho_a_yp = tr_rho_a_y_->pointer();
       double * tr_rho_b_yp = tr_rho_b_y_->pointer();
 
       double * tr_rho_a_zp = tr_rho_a_z_->pointer();
       double * tr_rho_b_zp = tr_rho_b_z_->pointer();
       
       double * tr_sigma_aap = tr_sigma_aa_->pointer();
       double * tr_sigma_abp = tr_sigma_ab_->pointer();
       double * tr_sigma_bbp = tr_sigma_bb_->pointer();

       double * rho_a_xp = rho_a_x_->pointer();
       double * rho_b_xp = rho_b_x_->pointer();

       double * rho_a_yp = rho_a_y_->pointer();
       double * rho_b_yp = rho_b_y_->pointer();
 
       double * rho_a_zp = rho_a_z_->pointer();
       double * rho_b_zp = rho_b_z_->pointer();

       double * pi_xp = pi_x_->pointer();
       double * pi_yp = pi_y_->pointer();
       double * pi_zp = pi_z_->pointer();

       for (int p = 0; p < phi_points_; p++) {

           double rho_x = rho_a_xp[p] + rho_b_xp[p];
           double rho_y = rho_a_yp[p] + rho_b_yp[p];
           double rho_z = rho_a_zp[p] + rho_b_zp[p];

           double rho = rho_p[p];
           double pi = pi_p[p];
           double DelR = R_p[p] - R1;

           double zeta = 0.0;
           double R = 0.0;
 
           // if ( !(rho < tol) && !(pi < tol) ) {
	   if ( !(rho < tol) && !(pi < 0.0) ) {
           
               R = R_p[p];
               // R = tanh(R);
           
               if ( ((1.0 - R) > tol) && ( R < R0 ) ) {
           
                  zeta = sqrt(1.0 - R);
           
                  tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0) + (R * rho_x) / (2.0*zeta) - pi_xp[p] / (rho*zeta);
                  tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0) - (R * rho_x) / (2.0*zeta) + pi_xp[p] / (rho*zeta);
                  
                  tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0) + (R * rho_y) / (2.0*zeta) - pi_yp[p] / (rho*zeta);
                  tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0) - (R * rho_y) / (2.0*zeta) + pi_yp[p] / (rho*zeta);

                  tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0) + (R * rho_z) / (2.0*zeta) - pi_zp[p] / (rho*zeta);
                  tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0) - (R * rho_z) / (2.0*zeta) + pi_zp[p] / (rho*zeta);

               }else if( !(R < R0) && !(R > R1) ) {
           
                       zeta = A * pow(DelR, 5.0) + B * pow(DelR, 4.0) + C * pow(DelR, 3.0);

                       tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0) 
                                       + (A * pow(DelR, 4.0)) * ( (10.0 * pi_xp[p] / rho) - (5.0 * R * rho_x) )
                                       + (B * pow(DelR, 3.0)) * ( (8.0  * pi_xp[p] / rho) - (4.0 * R * rho_x) )
                                       + (C * pow(DelR, 2.0)) * ( (6.0  * pi_xp[p] / rho) - (3.0 * R * rho_x) );

                       tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0) 
                                       + (A * pow(DelR, 4.0)) * (-(10.0 * pi_xp[p] / rho) + (5.0 * R * rho_x) )
                                       + (B * pow(DelR, 3.0)) * (-(8.0  * pi_xp[p] / rho) + (4.0 * R * rho_x) )
                                       + (C * pow(DelR, 2.0)) * (-(6.0  * pi_xp[p] / rho) + (3.0 * R * rho_x) );
           
                       tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0) 
                                       + (A * pow(DelR, 4.0)) * ( (10.0 * pi_yp[p] / rho) - (5.0 * R * rho_y) )
                                       + (B * pow(DelR, 3.0)) * ( (8.0  * pi_yp[p] / rho) - (4.0 * R * rho_y) )
                                       + (C * pow(DelR, 2.0)) * ( (6.0  * pi_yp[p] / rho) - (3.0 * R * rho_y) );

                       tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0) 
                                       + (A * pow(DelR, 4.0)) * (-(10.0 * pi_yp[p] / rho) + (5.0 * R * rho_y) )
                                       + (B * pow(DelR, 3.0)) * (-(8.0  * pi_yp[p] / rho) + (4.0 * R * rho_y) )
                                       + (C * pow(DelR, 2.0)) * (-(6.0  * pi_yp[p] / rho) + (3.0 * R * rho_y) );

                       tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0) 
                                       + (A * pow(DelR, 4.0)) * ( (10.0 * pi_zp[p] / rho) - (5.0 * R * rho_z) )
                                       + (B * pow(DelR, 3.0)) * ( (8.0  * pi_zp[p] / rho) - (4.0 * R * rho_z) )
                                       + (C * pow(DelR, 2.0)) * ( (6.0  * pi_zp[p] / rho) - (3.0 * R * rho_z) );

                       tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0) 
                                       + (A * pow(DelR, 4.0)) * (-(10.0 * pi_zp[p] / rho) + (5.0 * R * rho_z) )
                                       + (B * pow(DelR, 3.0)) * (-(8.0  * pi_zp[p] / rho) + (4.0 * R * rho_z) )
                                       + (C * pow(DelR, 2.0)) * (-(6.0  * pi_zp[p] / rho) + (3.0 * R * rho_z) );
               }else if( R > R1 ) {
           
                       zeta = 0.0;

                       tr_rho_a_xp[p] = (1.0 + zeta) * (rho_x/2.0);
                       tr_rho_b_xp[p] = (1.0 - zeta) * (rho_x/2.0);

                       tr_rho_a_yp[p] = (1.0 + zeta) * (rho_y/2.0);
                       tr_rho_b_yp[p] = (1.0 - zeta) * (rho_y/2.0);

                       tr_rho_a_zp[p] = (1.0 + zeta) * (rho_z/2.0);
                       tr_rho_b_zp[p] = (1.0 - zeta) * (rho_z/2.0);
               }
           
           }else{
           
               tr_rho_a_xp[p] = 0.0;
               tr_rho_b_xp[p] = 0.0;
               
               tr_rho_a_yp[p] = 0.0;
               tr_rho_b_yp[p] = 0.0;
               
               tr_rho_a_zp[p] = 0.0;
               tr_rho_b_zp[p] = 0.0;
           }
           // outfile->Printf("\n     Fully translated density gradient (x-component) =      %12.15lf\n",tr_rho_a_xp[p]);

           tr_sigma_aap[p] = (tr_rho_a_xp[p] * tr_rho_a_xp[p]) + (tr_rho_a_yp[p] * tr_rho_a_yp[p]) + (tr_rho_a_zp[p] * tr_rho_a_zp[p]);  
           tr_sigma_abp[p] = (tr_rho_a_xp[p] * tr_rho_b_xp[p]) + (tr_rho_a_yp[p] * tr_rho_b_yp[p]) + (tr_rho_a_zp[p] * tr_rho_b_zp[p]);  
           tr_sigma_bbp[p] = (tr_rho_b_xp[p] * tr_rho_b_xp[p]) + (tr_rho_b_yp[p] * tr_rho_b_yp[p]) + (tr_rho_b_zp[p] * tr_rho_b_zp[p]);  

        }

    }

}

void MCPDFTSolver::polyradical_analysis() {

        int mult = molecule_->multiplicity();


        if (mult != 1) {
           outfile->Printf("\n");
           outfile->Printf("    >>> For open-shell systems, the nat_orbs option\n");
           outfile->Printf("        in v2RDM-CASSCF block should be false. <<<\n");
        } 

        GetGridInfo();

        std::shared_ptr<PSIO> psio (new PSIO());

        // psio->set_pid("18332");

        if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
        if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

        // D1a

        psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

        long int na;
        psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

        opdm_a_ = (opdm *)malloc(na * sizeof(opdm));
        psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a_,na * sizeof(opdm));
        psio->close(PSIF_V2RDM_D1A,1);

        for (int n = 0; n < na; n++) {

            int i = opdm_a_[n].i;
            int j = opdm_a_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            if ( hi != hj ) {
                throw PsiException("error: something is wrong with the symmetry of the alpha OPDM",__FILE__,__LINE__);
            }

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            Da_->pointer(hi)[ii][jj] = opdm_a_[n].val;

        }

        // D1b

        psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

        long int nb;
        psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

        opdm_b_ = (opdm *)malloc(nb * sizeof(opdm));
        psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b_,nb * sizeof(opdm));
        psio->close(PSIF_V2RDM_D1B,1);

        for (int n = 0; n < nb; n++) {

            int i = opdm_b_[n].i;
            int j = opdm_b_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            if ( hi != hj ) {
                throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
            }

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            Db_->pointer(hi)[ii][jj] = opdm_b_[n].val;

        }
        std::shared_ptr<Matrix> Dtot;
        Dtot = std::shared_ptr<Matrix>(new Matrix(Da_));
        Dtot->add(Db_);
        // Dtot->print();

        std::shared_ptr<Matrix> Dmat (new Matrix(Dtot));
        std::shared_ptr<Matrix> Umat (new Matrix(Dtot));
        std::shared_ptr<Matrix> Nmat (new Matrix(Dtot));
        Umat->zero();
        Dmat->zero();
        Nmat->zero();
        for (int n = 0; n < na; n++) {

            int i = opdm_a_[n].i;
            int j = opdm_a_[n].j;

            int hi = symmetry_[i];
            int hj = symmetry_[j];

            if ( hi != hj ) {
                throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
            }

            int ii = i - pitzer_offset_[hi];
            int jj = j - pitzer_offset_[hj];

            double popval = Dtot->pointer(hi)[ii][jj];

            Umat->pointer(hi)[ii][jj] = 1.0 - abs(1.0 - popval);
            Dmat->pointer(hi)[ii][jj] = popval * (2.0 - popval);
            Nmat->pointer(hi)[ii][jj] = popval * popval * (2.0 - popval) * (2.0 - popval);
            
        }
        // Umat->back_transform(Ca_);
        // Umat->print(); 

        rho_    = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_a_  = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_  = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        Ur_     = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        Dr_     = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        Nr_     = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        double tempr  = 0.0;
        double tempra = 0.0;
        double temprb = 0.0;
        double tempu  = 0.0;
        double tempd  = 0.0;
        double tempn  = 0.0;
        for (int p = 0; p < phi_points_; p++) {
            double dumr  = 0.0;
            double dumra = 0.0;
            double dumrb = 0.0;
            double dumu  = 0.0;
            double dumd  = 0.0;
            double dumn  = 0.0;
            for (int n = 0; n < na; n++) {

                int i = opdm_a_[n].i;
                int j = opdm_a_[n].j;

                int hi = symmetry_[i];
                int hj = symmetry_[j];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];

                dumra += super_phi_->pointer(hi)[p][ii] *
                         super_phi_->pointer(hj)[p][jj] * opdm_a_[n].val;

                dumu += super_phi_->pointer(hi)[p][ii] *
                        super_phi_->pointer(hj)[p][jj] * Umat->pointer(hi)[ii][jj];

                dumd += super_phi_->pointer(hi)[p][ii] *
                        super_phi_->pointer(hj)[p][jj] * Dmat->pointer(hi)[ii][jj];

                dumn += super_phi_->pointer(hi)[p][ii] *
                        super_phi_->pointer(hj)[p][jj] * Nmat->pointer(hi)[ii][jj];
            }
            Ur_ ->pointer()[p] = dumu;
            Dr_ ->pointer()[p] = dumd;
            Nr_ ->pointer()[p] = dumn;

            for (int n = 0; n < nb; n++) {

                int i = opdm_b_[n].i;
                int j = opdm_a_[n].j;

                int hi = symmetry_[i];
                int hj = symmetry_[j];

                int ii = i - pitzer_offset_[hi];
                int jj = j - pitzer_offset_[hj];

                dumrb += super_phi_->pointer(hi)[p][ii] *
                         super_phi_->pointer(hj)[p][jj] * opdm_b_[n].val;
            }
            rho_  ->pointer()[p] += dumra + dumrb;
            rho_a_->pointer()[p] += dumra;
            rho_b_->pointer()[p] += dumrb;

            tempr  += grid_w_->pointer()[p] * (dumra + dumrb);
            tempra += grid_w_->pointer()[p] * dumra;
            temprb += grid_w_->pointer()[p] * dumrb;
            tempu  += grid_w_->pointer()[p] * dumu;
            tempd  += grid_w_->pointer()[p] * dumd;
            tempn  += grid_w_->pointer()[p] * dumn;
        }
        // Ur_->print();
        // Dr_->print();

        outfile->Printf("\n\n");
        outfile->Printf("    Trace[ rho(r)  ]:     %-20.12lf\n",tempr);
        outfile->Printf("    Trace[ rhoa(r) ]:     %-20.12lf\n",tempra);
        outfile->Printf("    Trace[ rhob(r) ]:     %-20.12lf\n",temprb);
        outfile->Printf("    Trace[ U(r)    ]:     %-20.12lf\n",tempu);
        outfile->Printf("    Trace[ D(r)    ]:     %-20.12lf\n",tempd);
        outfile->Printf("    Trace[ N(r)    ]:     %-20.12lf\n",tempn);

        double * w_p    = grid_w_->pointer(); 
        double * x_p    = grid_x_->pointer(); 
        double * y_p    = grid_y_->pointer(); 
        double * z_p    = grid_z_->pointer();
        double * rho_p  = rho_   ->pointer(); 
        double * rhoa_p = rho_a_ ->pointer(); 
        double * rhob_p = rho_b_ ->pointer(); 
        double * Ur_p   = Ur_    ->pointer(); 
        double * Dr_p   = Dr_    ->pointer(); 
        double * Nr_p   = Nr_    ->pointer(); 

        FILE *pfile;
        pfile = fopen("mygrids.txt","w");
        if (mult == 1) {
           std::fprintf(pfile,"w                    x                   y                  z                  rho(r)                  U(r)                  D(r)                    N(r)\n"); 
           for (int p = 0; p < phi_points_; p++) {
                 std::fprintf(pfile,"%-16.12lf             %-16.12lf            %-16.12lf            %-16.12lf             %-16.12lf             %-16.12lf             %-16.12lf              %-16.12lf\n"
                 ,w_p[p],x_p[p],y_p[p],z_p[p],rho_p[p],Ur_p[p],Dr_p[p],Nr_p[p]);
           }
        }else{
             std::fprintf(pfile,"w                    x                   y                  z                  rho(r)                 rho_s(r)                   U(r)                  D(r)                    N(r)\n"); 
             for (int p = 0; p < phi_points_; p++) {
                   std::fprintf(pfile,"%-16.12lf             %-16.12lf            %-16.12lf            %-16.12lf            %-16.12lf             %-16.12lf             %-16.12lf             %-16.12lf              %-16.12lf\n"
                   ,w_p[p],x_p[p],y_p[p],z_p[p],rho_p[p],(rhoa_p[p]-rhob_p[p]),Ur_p[p],Dr_p[p],Nr_p[p]);
             }
        }
        fclose(pfile);

        free(opdm_a_);
        free(opdm_b_);

        if (options_.get_bool("DO_UNPAIRED_MULLIKEN")) {
           Da_->zero();
           Db_->zero();
           if (options_.get_str("UNPAIRED_MULLIKEN_TYPE") == "N") { 
              Da_->add(Nmat);
           }else if (options_.get_str("UNPAIRED_MULLIKEN_TYPE") == "U") {
              Da_->add(Umat);
           }else if (options_.get_str("UNPAIRED_MULLIKEN_TYPE") == "D") {
              Da_->add(Dmat);
           }
        }
}

}} // end of namespaces

