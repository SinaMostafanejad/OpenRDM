/*
 *@BEGIN LICENSE
 *
 * MCPDFT, a plugin to:
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

#include <psi4/psi4-dec.h>
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libpsi4util/PsiOutStream.h>

#include "mcpdft_solver.h"

using namespace psi;



namespace psi{namespace mcpdft{

void MCPDFTSolver::ReadAllIntegrals(iwlbuf *Buf) {

    // get symmetry information

    symmetry_ = (int*)malloc(nmo_*sizeof(int));
    int count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < nmopi_[h]; i++) {
            symmetry_[count++] = h;
        }
    }

    // geminals, by symmetry

    for (int h = 0; h < nirrep_; h++) {
        std::vector < std::pair<int,int> > mygems;
        for (int i = 0; i < nmo_; i++) {
            for (int j = i; j < nmo_; j++) {
                int sym = symmetry_[i] ^ symmetry_[j];
                if (h==sym) {
                    mygems.push_back(std::make_pair(j,i));
                }

            }
        }
        gems_.push_back(mygems);
    }

    // geminal / orbital maps

    bas_  = (int***)malloc(nirrep_*sizeof(int**));
    ibas_ = (int***)malloc(nirrep_*sizeof(int**));

    for (int h = 0; h < nirrep_; h++) {

        ibas_[h] = (int**)malloc(nmo_*sizeof(int*));
        bas_[h]  = (int**)malloc(nmo_*nmo_*sizeof(int*));

        for (int i = 0; i < nmo_; i++) {
            ibas_[h][i] = (int*)malloc(nmo_*sizeof(int));
            for (int j = 0; j < nmo_; j++) {
                ibas_[h][i][j] = -999;
            }
        }

        for (int i = 0; i < nmo_*nmo_; i++) {
            bas_[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_[h][i][j] = -999;
            }
        }

        for (int n = 0; n < gems_[h].size(); n++) {
            int i = gems_[h][n].first;
            int j = gems_[h][n].second;

            ibas_[h][i][j] = n;
            ibas_[h][j][i] = n;
            bas_[h][n][0]  = i;
            bas_[h][n][1]  = j;
        }

    }

    // allocate memory for erfc integrals

    long int dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        dim += (long int)gems_[h].size() * ( (long int)gems_[h].size() + 1 ) / 2;
    }
    erfc_tei_ = (double*)malloc(dim*sizeof(double));
    memset((void*)erfc_tei_,'\0',dim*sizeof(double));

    // ok, i think we're ready to read the integrals from disk
    
    unsigned long int lastbuf;
    Label *lblptr;
    Value *valptr;
  
    unsigned long int idx, p, q, r, s, pq, rs, pqrs;
  
    lblptr = Buf->labels;
    valptr = Buf->values;
    lastbuf = Buf->lastbuf;
  
    outfile->Printf("\n");
    outfile->Printf("        Read integrals......");
    /**
      * first buffer (read in when Buf was initialized)
      */
    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
        p = (unsigned long int) lblptr[idx++];
        q = (unsigned long int) lblptr[idx++];
        r = (unsigned long int) lblptr[idx++];
        s = (unsigned long int) lblptr[idx++];
  
        double val = (double)valptr[Buf->idx];
  
        // none of this will work with frozen virtuals ...
  
        int hp = symmetry_[p];
        int hq = symmetry_[q];
        int hpq = hp ^ hq;
  
        pq = ibas_[hpq][p][q];
        rs = ibas_[hpq][r][s];
  
        long int offset = 0;
        for (int h = 0; h < hpq; h++) {
            offset += (long int)gems_[h].size() * ( (long int)gems_[h].size() + 1 ) / 2;
        }
        erfc_tei_[offset + INDEX(pq,rs)] = val;
  
    }
    /**
      * now do the same for the rest of the buffers
      */
    while(!lastbuf){
        iwl_buf_fetch(Buf);
        lastbuf = Buf->lastbuf;
        for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
  
            p = (unsigned long int) lblptr[idx++];
            q = (unsigned long int) lblptr[idx++];
            r = (unsigned long int) lblptr[idx++];
            s = (unsigned long int) lblptr[idx++];
  
            double val = (double)valptr[Buf->idx];
  
            int hp = symmetry_[p];
            int hq = symmetry_[q];
            int hpq = hp ^ hq;
  
            pq = ibas_[hpq][p][q];
            rs = ibas_[hpq][r][s];
  
            long int offset = 0;
            for (int h = 0; h < hpq; h++) {
                offset += (long int)gems_[h].size() * ( (long int)gems_[h].size() + 1 ) / 2;
            }
            erfc_tei_[offset + INDEX(pq,rs)] = val;
        }
  
    }
    outfile->Printf("done.\n\n");
}


void MCPDFTSolver::ReadRangeSeparatedIntegrals(){
    struct iwlbuf Buf;
    iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
    ReadAllIntegrals(&Buf);
    iwl_buf_close(&Buf,1);
}


}}
