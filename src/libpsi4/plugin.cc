/*
 * @BEGIN LICENSE
 *
 * mcpdft by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"
#include <psi4/libpsi4util/process.h>
#include "mcpdft_solver.h"

namespace psi{ namespace mcpdft {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "MCPDFT"|| options.read_globals()) {
        /*- MCPDFT type -*/
        options.add_str("MCPDFT_METHOD", "MCPDFT", "MCPDFT 1H_MCPDFT 1DH_MCPDFT RS_MCPDFT RS1H_MCPDFT RS1DH_MCPDFT");
        /*- The range-separation parameter -*/
        options.add_double("MCPDFT_OMEGA", 0.0);
        /*- Coupling parameter Lambda for hybrid MCPDFT functionals -*/
        options.add_double("MCPDFT_LAMBDA", 0.00);
        /*- Reference must be UKS -*/
        options.add_str("REFERENCE", "UKS");
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- MCPDFT functional -*/
        options.add_str("MCPDFT_FUNCTIONAL", "SVWN", "SVWN PBE REVPBE BOP BLYP WPBE LRC_WPBE");
        /*- type of density and density gradient translation:
        REGULAR = The gradients of on-top density are not considered in the polarization factor zeta
        FULL = The gradients of on-top density is included in the polarization factor zeta       -*/
        options.add_str("MCPDFT_TRANSLATION_TYPE", "REGULAR", "REGULAR FULL");
        /*- JK object type can be DF or PK -*/
        options.add_str("MCPDFT_TYPE", "DF", "DF PK");
        /*- reference type -*/
        options.add_str("MCPDFT_REFERENCE", "V2RDM", "V2RDM CI");
        /*- Construct effectively unpaired electron matrices, D(r) and U(r)? -*/
        options.add_bool("POLYRADICAL_ANALYSIS", false);
        /*- Type of the effectively unpaired electron population analysis -*/
        options.add_bool("DO_UNPAIRED_MULLIKEN", false);
        /*- Type of the effectively unpaired electron population analysis -*/
        options.add_str("UNPAIRED_MULLIKEN_TYPE","N", "N U D");

    }

    return true;
}

extern "C" PSI_API
SharedWavefunction mcpdft(SharedWavefunction ref_wfn, Options& options)
{

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

    //std::shared_ptr<MCPDFTSolver> dft (new MCPDFTSolver(ref_wfn,options));
    //double energy{0.0}; 
    //if (options.get_bool("POLYRADICAL_ANALYSIS"))
    //   dft->polyradical_analysis();
    //else
    //   energy = dft->compute_energy();

    //Process::environment.globals["CURRENT ENERGY"] = energy;
    //dft = NULL;
    // TODO: return mcpdft wave function instead of reference
    return ref_wfn;
}

}} // End namespaces

