import numpy as np
import scipy as sp
from functools import reduce
import time
import os
import sys
import tempfile
import h5py
import pyscf
from pyscf import gto, scf

from pyscf import mcscf
from pyscf.mcscf import addons

from pyscf.dmrgscf import dmrgci

from pyscf import dft
from pyscf.dft import numint

class MCPDFT:

   def __init__(self, mol, mc, ref_method=None):
      self.ref_method = ref_method
      self.cas = mc
      if mol == None:
         self.mol = gto.Mole()
      else:
         self.mol = mol
      # Check settings for DMRG
      if (self.ref_method == 'DMRG'):
         try:
             from pyscf.dmrgscf import settings
         except ImportError:
             settings = lambda: None
             settings.BLOCKEXE = getattr(__config__, 'dmrgscf_BLOCKEXE', None)
             settings.BLOCKEXE_COMPRESS_NEVPT = \
                     getattr(__config__, 'dmrgscf_BLOCKEXE_COMPRESS_NEVPT', None)
             settings.BLOCKSCRATCHDIR = getattr(__config__, 'dmrgscf_BLOCKSCRATCHDIR', None)
             settings.BLOCKRUNTIMEDIR = getattr(__config__, 'dmrgscf_BLOCKRUNTIMEDIR', None)
             settings.MPIPREFIX = getattr(__config__, 'dmrgscf_MPIPREFIX', None)
             settings.BLOCKVERSION = getattr(__config__, 'dmrgscf_BLOCKVERSION', None)
             if (settings.BLOCKEXE is None or settings.BLOCKSCRATCHDIR is None):
                 import sys
                 sys.stderr.write('settings.py not found for module dmrgci.  Please create %s\n'
                                  % os.path.join(os.path.dirname(__file__), 'settings.py'))
                 raise ImportError('settings.py not found')
 
         #self.cas.fcisolver = dmrgci.DMRGCI(mol, maxM=2000, tol=1.e-8)
         #self.cas.callback = self.cas.fcisolver.restart_scheduler_()
         #if self.cas.chkfile == self.cas._scf._chkfile.name:
         #   # Do not delete chkfile after mcscf
         #   self.cas.chkfile = tempfile.mktemp(dir=settings.BLOCKSCRATCHDIR)
         #   if not os.path.exists(settings.BLOCKSCRATCHDIR):
         #      os.makedirs(settings.BLOCKSCRATCHDIR)

      self.ci      = mc.ci
      self.C_mo    = mc.mo_coeff
      self.nao     = mc.mo_coeff.shape[0]
      self.nmo     = mc.mo_coeff.shape[1]
      self.nfrz    = mc.frozen
      self.ncore   = mc.ncore
      self.ncas    = mc.ncas
      self.nelecas = mc.nelecas

      self.nocc    = self.ncore + self.ncas
      self.virt    = self.nmo - self.nocc
      self.amo     = self.nmo - self.ncore - (0 if self.nfrz == None else self.nfrz)

      self.write_hdf5 = False

   def _print_active_space(self):
       print("\n")
       print("----------------------------------------------")
       print("   Number of AOs:                     %s" % self.nao)
       print("   Number of MOs:                     %s" % self.nmo)
       print("   Number of frozen orbitals:         %s" % ('0' if self.nfrz == None else self.nfrz) )
       print("   Number of core orbitals:           %s" % self.ncore)
       print("   Number of active orbitals:         %s" % self.ncas)
       print("   Number of active alpha electrons:  %s" % self.nelecas[0])
       print("   Number of active beta electrons:   %s" % self.nelecas[1])
       print("----------------------------------------------")
       print("\n")

   def make_active_rdm1s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      ncas    = cas.ncas
      nelecas = cas.nelecas 
      if self.ref_method == 'MCSCF':
         if ci  is None: ci  = self.ci
         casdm1a, casdm1b = cas.fcisolver.make_rdm1s(ci, ncas, nelecas)
      else: # ref_method == 'DMRG'
         # first argument takes 0 for ground and 1 for excited state
         casdm1a, casdm1b = cas.fcisolver.make_rdm1s(0, ncas, nelecas)
      return casdm1a, casdm1b

   def make_full_rdm1s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      nmo     = cas.mo_coeff.shape[1]
      ncore   = cas.ncore
      ncas    = cas.ncas
      nocc    = ncore + ncas
      nelecas = cas.nelecas 
      # building core part
      dm1a = np.zeros((nmo,nmo))
      dm1b = np.zeros((nmo,nmo))
      idx  = np.arange(ncore)
      dm1a[idx,idx] = dm1b[idx,idx] = 1.0
      # building active part
      casdm1a, casdm1b = self.make_active_rdm1s(cas=cas, ci=ci)
      dm1a[ncore:nocc,ncore:nocc] = casdm1a
      dm1b[ncore:nocc,ncore:nocc] = casdm1b
      return dm1a, dm1b

   def make_active_rdm2s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      casdm1s, casdm2s = cas.fcisolver.make_rdm12s(ci, cas.ncas, cas.nelecas, reorder=True)
      # aa, ab, bb
      return casdm2s[0], casdm2s[1], casdm2s[2]

   def make_full_rdm2s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      nmo     = cas.mo_coeff.shape[1]
      ncore   = cas.ncore
      ncas    = cas.ncas
      nocc    = ncore + ncas
      nelecas = cas.nelecas 
      # TODO: check to see if it is necessary to change nmo -> nocc
      dm2aa = np.zeros((nmo,nmo,nmo,nmo))
      dm2bb = np.zeros((nmo,nmo,nmo,nmo))
      dm2ab = np.zeros((nmo,nmo,nmo,nmo))
      casdm1a, casdm1b = self.make_active_rdm1s(cas=cas, ci=ci)
      casdm2aa, casdm2ab, casdm2bb = self.make_active_rdm2s()
      # Be aware that Chemist's notation should be adopted!
      #----------------
      # active-active part
      #----------------
      dm2aa[ncore:nocc,ncore:nocc,ncore:nocc,ncore:nocc] = casdm2aa
      dm2bb[ncore:nocc,ncore:nocc,ncore:nocc,ncore:nocc] = casdm2bb
      dm2ab[ncore:nocc,ncore:nocc,ncore:nocc,ncore:nocc] = casdm2ab
      for i in range(ncore):
         for j in range(ncore):
            #----------------
            # core-core part
            #----------------
            dm2aa[i,i,j,j] = \
            dm2bb[i,i,j,j] = \
            dm2ab[i,i,j,j] = 1.0

            dm2aa[i,j,j,i] = \
            dm2bb[i,j,j,i] = -1.0
            #----------------
            # core-active part
            #----------------
            ## aa block
            dm2aa[i,i,ncore:nocc,ncore:nocc] = \
            dm2aa[ncore:nocc,ncore:nocc,i,i] = casdm1a

            dm2aa[i,ncore:nocc,ncore:nocc,i] = \
            dm2aa[ncore:nocc,i,i,ncore:nocc] = -casdm1a

            ## bb block
            dm2bb[i,i,ncore:nocc,ncore:nocc] = \
            dm2bb[ncore:nocc,ncore:nocc,i,i] = casdm1b

            dm2bb[i,ncore:nocc,ncore:nocc,i] = \
            dm2bb[ncore:nocc,i,i,ncore:nocc] = -casdm1b

            ## ab block
            dm2ab[i,i,ncore:nocc,ncore:nocc] = casdm1b
            dm2ab[ncore:nocc,ncore:nocc,i,i] = casdm1a

      return dm2aa, dm2ab, dm2bb

   def make_active_rdm12(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      ncas    = cas.ncas
      nelecas = cas.nelecas
      # This class member can distinguish between DMRGSCF and MCSCF cas objects
      casdm1a, casdm1b = self.make_active_rdm1s(cas, ci)
      casdm1 = casdm1a + casdm1b
      if self.ref_method == 'MCSCF':
         casdm2aa, casdm2ab, casdm2bb = self.make_active_rdm2s(cas, ci)
         casdm2 = casdm2aa + casdm2ab + casdm2ab.transpose(2,3,0,1) + casdm2bb
      else: # ref_method == 'DMRG'
         casdm2 = cas.fcisolver.make_rdm12(0, ncas, nelecas)[1]
      return casdm1, casdm2

   def make_full_rdm12(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      nmo     = cas.mo_coeff.shape[1]
      ncore   = cas.ncore
      ncas    = cas.ncas
      nocc    = ncore + ncas
      nelecas = cas.nelecas 
      dm1 = np.zeros((nmo,nmo))
      dm2 = np.zeros((nmo,nmo,nmo,nmo))
      if self.ref_method == 'MCSCF':
         if ci  is None: ci  = self.ci
         dm1a , dm1b = self.make_full_rdm1s(cas, ci)
         dm1 = dm1a + dm1b
         dm2aa, dm2ab, dm2bb = self.make_full_rdm2s(cas,ci)
         dm2 = dm2aa + dm2ab + dm2ab.transpose(2,3,0,1) + dm2bb
      else: # ref_method == 'DMRG'
         #----------------
         # Be aware that Chemist's notation should be adopted!
         #----------------
         #   1-RDM
         #----------------
         # core part
         #----------------
         idx = np.arange(ncore)
         dm1[idx,idx] = 2
         #----------------
         # active part
         #----------------
         casdm1, casdm2 = cas.fcisolver.make_rdm12(0, ncas, nelecas)
         dm1[ncore:nocc,ncore:nocc] = casdm1
         #----------------
         #   2-RDM
         #----------------
         # active-active part
         #----------------
         dm2[ncore:nocc,ncore:nocc,ncore:nocc,ncore:nocc] = casdm2
         for i in range(ncore):
            for j in range(ncore):
               #----------------
               # core-core part
               #----------------
               dm2[i,i,j,j] = 4.0
               dm2[i,j,j,i] = -2.0
               #----------------
               # core-active part
               #----------------
               dm2[i,i,ncore:nocc,ncore:nocc] = \
               dm2[ncore:nocc,ncore:nocc,i,i] = 2.0*casdm1

               dm2[i,ncore:nocc,ncore:nocc,i] = \
               dm2[ncore:nocc,i,i,ncore:nocc] = -casdm1

      return dm1, dm2

   def write_rdm1s_d2sp_coo(self, dm1a=None, dm1b=None, is_active=False, f=None):
      row = dm1a.shape[0]
      col = dm1a.shape[1]
      assert(row == col)

      if dm1a.all() == None:
         if is_active == False:
            dm1a = self.dm1a
         else:
            dm1a = self.casdm1a

      if dm1b.all() == None:
         if is_active == False:
            dm1b = self.dm1b
         else:
            dm1b = self.casdm1b

      if f == None: f = h5py.File('data.h5','w')
  
      # TODO: maybe that's a good idea to make it a class member variable and
      # give the user the choice to choose the tolerance cutoff value
      tol = 1.0e-20
      # storing nonzero elements in the upper triangular part (i <= j)
      # row-based packing in COO format
      val     = []
      row_idx = []
      col_idx = []
      nnz = 0  #number of non-zero elements in the upper triangle

      # alpha block
      for i in range(col):
         for j in range(i,col):
            dum = dm1a[i,j]
            if (abs(dum) > tol):
               val.append(dum)
               row_idx.append(i)
               col_idx.append(j)
               nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D1/FULL_D1a_MO/NNZ"] = nnz
         f["/SP_SYM_D1/FULL_D1a_MO/VAL"] = val
         f["/SP_SYM_D1/FULL_D1a_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/FULL_D1a_MO/COL_IDX"] = col_idx
      else:
         f["/SP_SYM_D1/ACT_D1a_MO/NNZ"] = nnz
         f["/SP_SYM_D1/ACT_D1a_MO/VAL"] = val
         f["/SP_SYM_D1/ACT_D1a_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/ACT_D1a_MO/COL_IDX"] = col_idx

      dum = 0.0
      val.clear()
      row_idx.clear()
      col_idx.clear()
      nnz = 0  #number of non-zero elements in the upper triangle

      # beta block
      for i in range(col):
         for j in range(i,col):
            dum = dm1b[i,j]
            if (abs(dum) > tol):
               val.append(dum)
               row_idx.append(i)
               col_idx.append(j)
               nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D1/FULL_D1b_MO/NNZ"] = nnz
         f["/SP_SYM_D1/FULL_D1b_MO/VAL"] = val
         f["/SP_SYM_D1/FULL_D1b_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/FULL_D1b_MO/COL_IDX"] = col_idx
      else:
         f["/SP_SYM_D1/ACT_D1b_MO/NNZ"] = nnz
         f["/SP_SYM_D1/ACT_D1b_MO/VAL"] = val
         f["/SP_SYM_D1/ACT_D1b_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/ACT_D1b_MO/COL_IDX"] = col_idx

   def write_rdm2s_d2sp_coo(self, dm2aa=None, dm2ab=None, dm2bb=None, is_active=False, f=None):
      dim1 = dm2aa.shape[0]
      dim2 = dm2aa.shape[1]
      dim3 = dm2aa.shape[2]
      dim4 = dm2aa.shape[3]

      assert(dim1 == dim2 == dim3 == dim4)
      if dm2aa.all() == None:
         if is_active == False:
            dm2aa = self.dm2aa
         else:
            dm2aa = self.casdm2aa

      if dm2ab.all() == None:
         if is_active == False:
            dm2ab = self.dm2ab
         else:
            dm2ab = self.casdm2ab

      if dm2bb.all() == None:
         if is_active == False:
            dm2bb = self.dm2bb
         else:
            dm2bb = self.casdm2bb

      if f == None: f = h5py.File('data.h5','w')
   
      # TODO: maybe that's a good idea to make it a class member variable and
      # give the user the choice to choose the tolerance cutoff value
      tol = 1.0e-20
      # storing nonzero elements in the upper triangular part (i <= j)
      # row-based packing in COO format
      val      = []
      dim1_idx = []
      dim2_idx = []
      dim3_idx = []
      dim4_idx = []
      nnz = 0  #number of non-zero elements in the upper triangle

      # alpha-alpha block
      for i in range(dim1):
         for j in range(i,dim2):
            for k in range(dim3):
               for l in range(k,dim4):
                  dum = dm2aa[i,j,k,l]
                  if (abs(dum) > tol):
                     val.append(dum)
                     dim1_idx.append(i)
                     dim2_idx.append(j)
                     dim3_idx.append(k)
                     dim4_idx.append(l)
                     nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D2/FULL_D2aa_MO/NNZ"] = nnz
         f["/SP_SYM_D2/FULL_D2aa_MO/VAL"] = val
         f["/SP_SYM_D2/FULL_D2aa_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/FULL_D2aa_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/FULL_D2aa_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/FULL_D2aa_MO/DIM4_IDX"] = dim4_idx
      else:
         f["/SP_SYM_D2/ACT_D2aa_MO/NNZ"] = nnz
         f["/SP_SYM_D2/ACT_D2aa_MO/VAL"] = val
         f["/SP_SYM_D2/ACT_D2aa_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/ACT_D2aa_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/ACT_D2aa_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/ACT_D2aa_MO/DIM4_IDX"] = dim4_idx

      dum = 0.0
      val.clear()
      dim1_idx.clear() 
      dim2_idx.clear()
      dim3_idx.clear()
      dim4_idx.clear()
      nnz = 0  #number of non-zero elements in the upper triangle

      for i in range(dim1):
         for j in range(dim2):
            for k in range(dim3):
               for l in range(dim4):
                  dum = dm2ab[i,j,k,l]
                  if (abs(dum) > tol):
                     val.append(dum)
                     dim1_idx.append(i)
                     dim2_idx.append(j)
                     dim3_idx.append(k)
                     dim4_idx.append(l)
                     nnz = nnz + 1
#                     print("val, dim1, dim2, dim3, dim4 = {}, {}, {}, {}, {}\n" .format(val,i,j,k,l))

#      # alpha-beta block
#      for i in range(dim1):
#         for j in range(i,dim2):
#            for k in range(dim3):
#               for l in range(k,dim4):
#                  dum = dm2ab[i,j,k,l]
#                  if (abs(dum) > tol):
#                     val.append(dum)
#                     dim1_idx.append(i)
#                     dim2_idx.append(j)
#                     dim3_idx.append(k)
#                     dim4_idx.append(l)
#                     nnz = nnz + 1
#                     print("val, dim1, dim2, dim3, dim4 = {}, {}, {}, {}, {}\n" .format(val,i,j,k,l))

      if (is_active == False):
         f["/SP_SYM_D2/FULL_D2ab_MO/NNZ"] = nnz
         f["/SP_SYM_D2/FULL_D2ab_MO/VAL"] = val
         f["/SP_SYM_D2/FULL_D2ab_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/FULL_D2ab_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/FULL_D2ab_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/FULL_D2ab_MO/DIM4_IDX"] = dim4_idx
      else:
         f["/SP_SYM_D2/ACT_D2ab_MO/NNZ"] = nnz
         f["/SP_SYM_D2/ACT_D2ab_MO/VAL"] = val
         f["/SP_SYM_D2/ACT_D2ab_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/ACT_D2ab_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/ACT_D2ab_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/ACT_D2ab_MO/DIM4_IDX"] = dim4_idx

      dum = 0.0
      val.clear()
      dim1_idx.clear() 
      dim2_idx.clear()
      dim3_idx.clear()
      dim4_idx.clear()
      nnz = 0  #number of non-zero elements in the upper triangle

      # beta-beta block
      for i in range(dim1):
         for j in range(i,dim2):
            for k in range(dim3):
               for l in range(k,dim4):
                  dum = dm2bb[i,j,k,l]
                  if (abs(dum) > tol):
                     val.append(dum)
                     dim1_idx.append(i)
                     dim2_idx.append(j)
                     dim3_idx.append(k)
                     dim4_idx.append(l)
                     nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D2/FULL_D2bb_MO/NNZ"] = nnz
         f["/SP_SYM_D2/FULL_D2bb_MO/VAL"] = val
         f["/SP_SYM_D2/FULL_D2bb_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/FULL_D2bb_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/FULL_D2bb_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/FULL_D2bb_MO/DIM4_IDX"] = dim4_idx
      else:
         f["/SP_SYM_D2/ACT_D2bb_MO/NNZ"] = nnz
         f["/SP_SYM_D2/ACT_D2bb_MO/VAL"] = val
         f["/SP_SYM_D2/ACT_D2bb_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/ACT_D2bb_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/ACT_D2bb_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/ACT_D2bb_MO/DIM4_IDX"] = dim4_idx

   def write_rdm12_d2sp_coo(self, dm1=None, dm2=None, is_active=False, f=None):
      row = dm1.shape[0]
      col = dm1.shape[1]

      dim1 = dm2.shape[0]
      dim2 = dm2.shape[1]
      dim3 = dm2.shape[2]
      dim4 = dm2.shape[3]

      assert(row == col)
      assert(dim1 == dim2 == dim3 == dim4)

      if dm1.all() == None:
         if is_active == False:
            dm1 = self.dm1
         else:
            dm1 = self.casdm1

      if dm2.all() == None:
         if is_active == False:
            dm2 = self.dm2
         else:
            dm2 = self.casdm2

      if f == None: f = h5py.File('data.h5','w')
  
      # TODO: maybe that's a good idea to make it a class member variable and
      # give the user the choice to choose the tolerance cutoff value
      tol = 1.0e-20
      # storing nonzero elements in the upper triangular part (i <= j)
      # row-based packing in COO format
      val     = []
      row_idx = []
      col_idx = []
      nnz = 0  #number of non-zero elements in the upper triangle

      # writing the non-zero elements of 1-RDM
      for i in range(col):
         for j in range(i,col):
            dum = dm1[i,j]
            if (abs(dum) > tol):
               val.append(dum)
               row_idx.append(i)
               col_idx.append(j)
               nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D1/FULL_D1_MO/NNZ"] = nnz
         f["/SP_SYM_D1/FULL_D1_MO/VAL"] = val
         f["/SP_SYM_D1/FULL_D1_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/FULL_D1_MO/COL_IDX"] = col_idx
      else:
         f["/SP_SYM_D1/ACT_D1_MO/NNZ"] = nnz
         f["/SP_SYM_D1/ACT_D1_MO/VAL"] = val
         f["/SP_SYM_D1/ACT_D1_MO/ROW_IDX"] = row_idx
         f["/SP_SYM_D1/ACT_D1_MO/COL_IDX"] = col_idx

      dum = 0.0
      val.clear()
      dim1_idx = []
      dim2_idx = []
      dim3_idx = []
      dim4_idx = []
      nnz = 0  #number of non-zero elements in the upper triangle

      # writing the non-zero elements of 2-RDM
      for i in range(dim1):
         for j in range(i,dim2):
            for k in range(dim3):
               for l in range(k,dim4):
                  dum = dm2[i,j,k,l]
                  if (abs(dum) > tol):
                     val.append(dum)
                     dim1_idx.append(i)
                     dim2_idx.append(j)
                     dim3_idx.append(k)
                     dim4_idx.append(l)
                     nnz = nnz + 1

      if (is_active == False):
         f["/SP_SYM_D2/FULL_D2_MO/NNZ"] = nnz
         f["/SP_SYM_D2/FULL_D2_MO/VAL"] = val
         f["/SP_SYM_D2/FULL_D2_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/FULL_D2_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/FULL_D2_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/FULL_D2_MO/DIM4_IDX"] = dim4_idx
      else:
         f["/SP_SYM_D2/ACT_D2_MO/NNZ"] = nnz
         f["/SP_SYM_D2/ACT_D2_MO/VAL"] = val
         f["/SP_SYM_D2/ACT_D2_MO/DIM1_IDX"] = dim1_idx
         f["/SP_SYM_D2/ACT_D2_MO/DIM2_IDX"] = dim2_idx
         f["/SP_SYM_D2/ACT_D2_MO/DIM3_IDX"] = dim3_idx
         f["/SP_SYM_D2/ACT_D2_MO/DIM4_IDX"] = dim4_idx

   def kernel(self):

      #--------------------------------------------- Active space info
      self._print_active_space()
      #--------------------------------------------- RDMs
      casdm1a, casdm1b = self.make_active_rdm1s()   # OK with both DMRG and MCSCF ref_methods
      #casdm1 = casdm1a + casdm1b
      casdm1 = self.make_active_rdm12()[0]          # OK with both DMRG and MCSCF ref_methods
      dm1a, dm1b = self.make_full_rdm1s()           # OK with both DMRG and MCSCF ref_methods
      #dm1 = dm1a + dm1b
      dm1 = self.make_full_rdm12()[0]               # OK with both DMRG and MCSCF ref_methods

      if (self.ref_method == 'MCSCF'):
         casdm2aa, casdm2ab, casdm2bb = self.make_active_rdm2s()    # Works only with ref_method == 'MCSCF'
         casdm2 = casdm2aa + casdm2ab + casdm2ab.transpose(2,3,0,1) + casdm2bb

         dm2aa, dm2ab, dm2bb = self.make_full_rdm2s()               # Works only with ref_method == 'MCSCF'
         dm2 = dm2aa + dm2ab + dm2ab.transpose(2,3,0,1) + dm2bb
      else: # ref_method == 'DMRG'
         casdm2 = self.make_active_rdm12()[1]       # OK with both DMRG and MCSCF ref_methods
         dm2  = self.make_full_rdm12()[1]           # OK with both DMRG and MCSCF ref_methods
      #print("\n D2aa:\n %s" % casdm2aa)
      #print("\n D2ab:\n %s" % casdm2ab)
      #print("\n D2bb:\n %s" % casdm2bb)
      #print("\n CAS D2:\n %s" % casdm2)
      #print("\n FULL D2:\n %s" % dm2)
      #print(casdm2 == dm2)
      #--------------------------------------------- nuclear repulsion 
      E_nn   = self.cas._scf.energy_nuc()
      #--------------------------------------------- core parts
      h      = self.cas._scf.get_hcore()
      h_mo   = ao2mo_transform(self.C_mo,h)
      E_core = hcore_energy(h_mo,dm1)
      #--------------------------------------------- J and K parts
      Ja = self.cas._scf.get_j (dm=ao2mo_transform(self.C_mo.T,dm1a))
      Jb = self.cas._scf.get_j (dm=ao2mo_transform(self.C_mo.T,dm1b))
      Ja_mo = ao2mo_transform(self.C_mo,Ja)
      Jb_mo = ao2mo_transform(self.C_mo,Jb)
      E_j   = Hartree_energy(Ja_mo, Jb_mo, dm1a, dm1b)
      #print(hcore_energies(h_mo, dm1a, dm1b))
      #---------------------------------------------

      print('=======================================================')
      print("   nuclear repulsion energy:     % .8f" % E_nn)
      print("   one-electron energy:          % .8f" % E_core)
      print("   classical Coulomb energy:     % .8f" % E_j)
      print('=======================================================')

      #--------------------------------------------- extracting grids and orbital values on them
      coords, weights, ao_values = get_grid_info(self.mol)
      mo_values = np.matmul(ao_values, self.C_mo)
      phi_ao, phi_ao_x, phi_ao_y, phi_ao_z = ao_values
      phi_mo, phi_mo_x, phi_mo_y, phi_mo_z = mo_values

      #--------------------------------------------- writing MCPDFT ingredients for OpenRDM into HDF5 file
      f = h5py.File("data.h5",'w')

      f["/N/N_AO"]      = self.nao
      f["/N/N_MO"]      = self.nmo
      f["/N/N_FRZ"]     = (0 if self.nfrz == None else self.nfrz)
      f["/N/N_COR"]     = self.ncore
      f["/N/N_CAS_ORB"] = self.ncas
      f["/N/N_CAS_ELE"] = self.nelecas
      f["/N/N_OCC_ACT"] = self.nocc
      f["/N/N_VIR_ACT"] = self.virt
      f["/N/N_MO_ACT"]  = self.amo
      f["/N/N_PTS"]     = coords.shape[0]

      f["/C"]        = self.C_mo

      f["/H/H_CORE_AO"] = h 
      f["/H/H_CORE_MO"] = h_mo

      f["/E/E_NUC"]     = E_nn 
      f["/E/E_CORE"]    = E_core 
      f["/E/E_HARTREE"] = E_j 

      f["/J/JA_AO"] = Ja
      f["/J/JB_AO"] = Jb
      f["/J/JA_MO"] = Ja_mo
      f["/J/JB_MO"] = Jb_mo

      f["/D/D1/FULL_D1A_MO"] = dm1a 
      f["/D/D1/FULL_D1B_MO"] = dm1b
      f["/D/D1/FULL_D1_MO"]  = dm1 

      f["/D/D1/ACT_D1A_MO"] = casdm1a 
      f["/D/D1/ACT_D1B_MO"] = casdm1b
      f["/D/D1/ACT_D1_MO"]  = casdm1 

      if self.ref_method == 'MCSCF':
         # turning these RDMs to physicist notation to be able to use the symmetry in sparsity
         dm2aa_r = dm2aa.transpose(0,2,1,3).reshape((self.nmo*self.nmo, self.nmo*self.nmo))
         #dm2aa_r = dm2aa.reshape((self.nmo*self.nmo, self.nmo*self.nmo))
         dm2ab_r = dm2ab.transpose(0,2,1,3).reshape((self.nmo*self.nmo, self.nmo*self.nmo))
         #dm2ab_r = dm2ab.reshape((self.nmo*self.nmo, self.nmo*self.nmo))
         dm2bb_r = dm2bb.transpose(0,2,1,3).reshape((self.nmo*self.nmo, self.nmo*self.nmo))
         #dm2bb_r = dm2bb.reshape((self.nmo*self.nmo, self.nmo*self.nmo))

         casdm2aa_r = casdm2aa.transpose(0,2,1,3).reshape((self.ncas*self.ncas, self.ncas*self.ncas))
         #casdm2aa_r = casdm2aa.reshape((self.ncas*self.ncas, self.ncas*self.ncas))
         casdm2ab_r = casdm2ab.transpose(0,2,1,3).reshape((self.ncas*self.ncas, self.ncas*self.ncas))
         #casdm2ab_r = casdm2ab.reshape((self.ncas*self.ncas, self.ncas*self.ncas))
         casdm2bb_r = casdm2bb.transpose(0,2,1,3).reshape((self.ncas*self.ncas, self.ncas*self.ncas))
         #casdm2bb_r = casdm2bb.reshape((self.ncas*self.ncas, self.ncas*self.ncas))

         f["/D/D2/FULL_D2AA_MO"] = dm2aa_r 
         f["/D/D2/FULL_D2AB_MO"] = dm2ab_r
         f["/D/D2/FULL_D2BB_MO"] = dm2bb_r

         f["/D/D2/ACT_D2AA_MO"] = casdm2aa_r 
         f["/D/D2/ACT_D2AB_MO"] = casdm2ab_r
         f["/D/D2/ACT_D2BB_MO"] = casdm2bb_r
      
      dm2_r = dm2.reshape((self.nmo*self.nmo, self.nmo*self.nmo))
      casdm2_r = casdm2.reshape((self.ncas*self.ncas, self.ncas*self.ncas))
      f["/D/D2/FULL_D2_MO"]  = dm2_r
      f["/D/D2/ACT_D2_MO"]   = casdm2_r

      #print(dm2ab)
      #print(dm2ab_r)

      f["/GRIDS/W"] = weights
      f["/GRIDS/X"] = coords[:,0]
      f["/GRIDS/Y"] = coords[:,1]
      f["/GRIDS/Z"] = coords[:,2]

      f["/PHI/PHI_AO"]   = phi_ao 
      f["/PHI/PHI_AO_X"] = phi_ao_x 
      f["/PHI/PHI_AO_Y"] = phi_ao_y 
      f["/PHI/PHI_AO_Z"] = phi_ao_z 

      f["/PHI/PHI_MO"]   = phi_mo 
      f["/PHI/PHI_MO_X"] = phi_mo_x 
      f["/PHI/PHI_MO_Y"] = phi_mo_y 
      f["/PHI/PHI_MO_Z"] = phi_mo_z 
      
      # writing spin blocks of full 1-RDMs into HDF5 file object f
      # (COO sparse format with matrix symmetry)
      self.write_rdm1s_d2sp_coo(dm1a,dm1b,is_active=False,f=f)
      # writing spin blocks of active-space 1-RDMs into HDF5 file object f
      # (COO sparse format with matrix symmetry)
      self.write_rdm1s_d2sp_coo(casdm1a,casdm1b,is_active=True,f=f)

      # writing spin blocks of active-space 2-RDMs into HDF5 file object f
      # (COO sparse format with matrix symmetry)
      if self.ref_method == 'MCSCF':
         self.write_rdm2s_d2sp_coo(casdm2aa,casdm2ab,casdm2bb,is_active=True,f=f)
         self.write_rdm2s_d2sp_coo(dm2aa,dm2ab,dm2bb,is_active=False,f=f)
      else: # ref_method == 'DMRG'
         # writing spin-free full 1- and 2-RDMs into HDF5 file object f
         # (COO sparse format with matrix symmetry)
         self.write_rdm12_d2sp_coo(dm1,dm2,is_active=False,f=f)
         # writing spin-free active-space 1- and 2-RDMs into HDF5 file object f
         # (COO sparse format with matrix symmetry)
         self.write_rdm12_d2sp_coo(casdm1,casdm2,is_active=True,f=f)

      return self

def ao2mo_transform(C_mo, mat):
   return reduce(np.dot,(C_mo.T,mat,C_mo))

def mo2ao_transform(C_mo, mat):
   return reduce(np.dot,(C_mo,mat,C_mo.T))

def hcore_energy(h, dm1):
   return reduce(np.tensordot,(h,dm1))

def hcore_energies(h, dm1a, dm1b):
   return reduce(np.tensordot,(h,dm1a)), \
          reduce(np.tensordot,(h,dm1b))
          
def Hartree_energy(Ja, Jb, dm1a, dm1b):
   Ej_aa = reduce(np.tensordot,(Ja,dm1a))
   Ej_ab = reduce(np.tensordot,(Ja,dm1b))
   Ej_ba = reduce(np.tensordot,(Jb,dm1a))
   Ej_bb = reduce(np.tensordot,(Jb,dm1b))
   return 0.5 * ( Ej_aa + Ej_ab + Ej_ba + Ej_bb )

def get_grid_info(mol=None):
   if mol == None:
      mol = gto.Mole()
   dft_obj = dft.RKS(mol)
   dft_obj.grids.atom_grid = (75, 302)
   dft_obj.grids.radi_method = dft.treutler
   dft_obj.grids.prune = dft.nwchem_prune
   #dft_obj.grids.radi_method = dft.mura_knowles
   #dft_obj.grids.becke_scheme = dft.original_becke
   #dft_obj.grids.level = 3
   dft_obj.kernel()
   # Orbital energies, Mulliken population etc.
   #dft_obj.analyze()

   # coords(n_points, 3): second dimension denotes x, y, z
   coords  = dft_obj.grids.coords
   weights = dft_obj.grids.weights
   #print(coords.shape)
   #print(weights.shape)

   # phi(4, n_points, nao): first dimension shows phi, phi_x, phi_y and phi_z, respectively
   ao_values = numint.eval_ao(mol, coords, deriv=1)  
   #print(phi)
   #print(phi_x)
   #print(ao_value.shape)

   # The first row of rho is electron density, the rest three rows are electron
   # density gradients which are needed for GGA functional
   #rho = numint.eval_rho(mol, ao_values, dm1, xctype='GGA')
   return coords, weights, ao_values




# ---------------junk----------------
#Te_Vne = np.tensordot (h, dm1)
#E_ja = np.dot(vj_mo, dm1_mo) / 2.0
#vj, vk = cas._scf.get_jk (dm=dm1s)
#vj = vj[0] + vj[1]
#print(np.allclose(dm1, dm1_alpha+dm1_beta))


# --------------RDM junk-------------
#dm1s = np.asarray ([dm1_alpha,dm1_beta])

#dm1 = self.cas.make_rdm1()
#dm1_alpha, dm1_beta = self.cas.make_rdm1s()
#print(dm1_alpha)
#print(dm1_beta)

#dm1_alpha, dm1_beta = pyscf.mcscf.addons.make_rdm1s(self.cas, mo_coeff=self.C_mo, ci=self.ci)
#print(dm1_alpha)
#print(dm1_beta)

#dm1_alpha_mo = self.ao2mo_transform(self.C_mo,dm1_alpha)
#dm1_beta_mo  = self.ao2mo_transform(self.C_mo,dm1_beta)
#print(dm1_alpha_mo)
#print(dm1_beta_mo)
#casdm1, casdm2 = self.cas.fcisolver.make_rdm12(self.ci, self.ncas, self.nelecas)
#dm1, dm2 = pyscf.mcscf.addons._make_rdm12_on_mo(casdm1, casdm2, self.ncore, self.ncas, self.nmo)
