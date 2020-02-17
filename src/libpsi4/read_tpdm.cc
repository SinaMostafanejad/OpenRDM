/*
 *@BEGIN LICENSE
 *
 * mcpdft, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/vector.h>

#include <psi4/libpsi4util/PsiOutStream.h>

#include <fstream>
#include <iostream>
#include <iomanip>

#include "mcpdft_solver.h"

using namespace psi;

namespace psi{namespace mcpdft{

void MCPDFTSolver::ReadTPDM() {

    std::shared_ptr<PSIO> psio (new PSIO());

    // TODO: should be added back when reading-in the density PSIOH files 
    // psio->set_pid ("18332");

    if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    
    //Ca_->print();

    psio_address addr_ab = PSIO_ZERO;

    // ab
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    long int nab;
    psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

    tpdm * d2 = (tpdm *)malloc(nab * sizeof(tpdm));
    memset((void*)d2,'\0',nab * sizeof(tpdm));

    psio->read_entry(PSIF_V2RDM_D2AB,"D2ab",(char*)d2,nab * sizeof(tpdm));

    psio->close(PSIF_V2RDM_D2AB,1);

    // build on-top pair density (already built of REFERENCE_TPDM = V2RDM)
    outfile->Printf("\n");
    outfile->Printf("    ==> Build Pi ...");
    if ( !is_low_memory_ ) {
        BuildPiFast(d2,nab);
    }else {
        BuildPiLowMemory(d2,nab);
    }
    outfile->Printf(" Done. <==\n");

/*
    // aa
    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);
    long int naa;
    psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));
    tpdm * d2aa = (tpdm *)malloc(naa * sizeof(tpdm));
    memset((void*)d2aa,'\0',naa * sizeof(tpdm));
    psio->read_entry(PSIF_V2RDM_D2AA,"D2aa",(char*)d2aa,naa * sizeof(tpdm));
    psio->close(PSIF_V2RDM_D2AA,1);

    // bb
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);
    long int nbb;
    psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));
    tpdm * d2bb = (tpdm *)malloc(nbb * sizeof(tpdm));
    memset((void*)d2bb,'\0',nbb * sizeof(tpdm));
    psio->read_entry(PSIF_V2RDM_D2BB,"D2bb",(char*)d2bb,nbb * sizeof(tpdm));
    psio->close(PSIF_V2RDM_D2BB,1);

    int myp = 0;
    double * x_p = grid_x_->pointer();
    double * y_p = grid_y_->pointer();
    double * z_p = grid_z_->pointer();
    for (int p = 0; p < phi_points_; p++) {
        if ( fabs(x_p[p] - 0.0) > 1e-6 ) continue;
        if ( fabs(y_p[p] - 0.0) > 1e-6 ) continue;
        if ( fabs(z_p[p] - 0.0) > 1e-6 ) continue;
        myp = p;
    }
 
    BuildExchangeCorrelationHole(myp,d2,nab,d2aa,naa,d2bb,nbb);
*/

    free(d2);

}

void MCPDFTSolver::ReadCITPDM(double* D, const char* fileName) {

    std::ifstream dataIn;

    dataIn.open(fileName);

    if (!dataIn) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    else {
         for (int i = 0; i < nmo_; i++)
             for (int j = 0; j < nmo_; j++)
                 for (int k = 0; k < nmo_; k++)
                     for (int l = 0; l < nmo_; l++) {

                         dataIn >> D[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                         if (D[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] < 1e-20)
                             D[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l] = 0.0;
                     }
    dataIn.close();
    }
}

void MCPDFTSolver::PrintTPDM(double* D) {

   for (int i = 0; i < nmo_; i++)
       for (int j = 0; j < nmo_; j++)
           for (int k = 0; k < nmo_; k++)
               for (int l = 0; l < nmo_; l++)
                   printf("%20.15lf\t", D[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);


   printf("\n\n");
}

}} //end namespaces


