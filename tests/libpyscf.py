import numpy as np
import scipy as sp
from functools import reduce
import time
import h5py
import pyscf
from pyscf import gto, scf, mcscf
from pyscf.mcscf import addons

class MCPDFT(mcscf.__class__):

   def __init__(self, mcscf):
      self.cas     = mcscf
      self.ci      = mcscf.ci
      self.C_mo    = mcscf.mo_coeff
      self.nao     = mcscf.mo_coeff.shape[0]
      self.nmo     = mcscf.mo_coeff.shape[1]
      self.nfrz    = mcscf.frozen
      self.ncore   = mcscf.ncore
      self.ncas    = mcscf.ncas
      self.nelecas = mcscf.nelecas

      self.nocc    = self.ncore + self.ncas
      self.virt    = self.nmo - self.nocc
      self.amo     = self.nmo - self.ncore - (0 if self.nfrz == None else self.nfrz)

   def print_active_space(self):
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


   def ao2mo_transform(self, C_mo, mat):
      return reduce(np.dot,(C_mo.T,mat,C_mo))
   
   def mo2ao_transform(self, C_mo, mat):
      return reduce(np.dot,(C_mo,mat,C_mo.T))
   
   def hcore_energy(self, h, dm1):
      return reduce(np.tensordot,(h,dm1))
   
   def hcore_energies(self, h, dm1a, dm1b):
      return reduce(np.tensordot,(h,dm1a)), \
             reduce(np.tensordot,(h,dm1b))
             
   def Hartree_energy(self, Ja, Jb, dm1a, dm1b):
      Ej_aa = reduce(np.tensordot,(Ja,dm1a))
      Ej_ab = reduce(np.tensordot,(Ja,dm1b))
      Ej_ba = reduce(np.tensordot,(Jb,dm1a))
      Ej_bb = reduce(np.tensordot,(Jb,dm1b))
      return 0.5 * ( Ej_aa + Ej_ab + Ej_ba + Ej_bb )

   def make_active_rdm1s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      casdm1a, casdm1b = cas.fcisolver.make_rdm1s(ci, cas.ncas, cas.nelecas)
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
      dm1a = dm1b = np.zeros((nmo,nmo))
      idx = np.arange(ncore)
      dm1a[idx,idx] = dm1b[idx,idx] = 1.0
      # building active part
      casdm1a, casdm1b = self.make_active_rdm1s(cas=cas, ci=ci)
      dm1a[ncore:nocc,ncore:nocc] = casdm1a
      dm1b[ncore:nocc,ncore:nocc] = casdm1b
      return dm1a, dm1b

   def make_active_rdm2s(self, cas=None, ci=None):
      if cas is None: cas = self.cas
      if ci  is None: ci  = self.ci
      nmo     = cas.mo_coeff.shape[1]
      ncore   = cas.ncore
      ncas    = cas.ncas
      nocc    = ncore + ncas
      nelecas = cas.nelecas 
      casdm1a, casdm1b = self.make_active_rdm1s(cas=cas, ci=ci)
      dm2aa = dm2bb, dm2ab = np.zeros((nmo,nmo,nmo,nmo))
      # Be aware that Chemist's notation should be adopted!
      for i in range(ncore):
         for j in range(ncore):
         dm2aa[i,i,ncore:nocc,ncore:nocc] = \
         dm2aa[ncore:nocc,ncore:nocc,i,i] = casdm1a
         dm2bb[i,ncore:nocc,ncore:nocc,i] = \
         dm2bb[ncore:nocc,i,i,ncore:nocc] = casdm1b
      return d2aa, d2bb, d2ab

   def kernel(self):

      #--------------------------------------------- Active space info
      self.print_active_space()
      #--------------------------------------------- RDMs
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
      #print(dm1)

      casdm1a, casdm1b = self.make_active_rdm1s()
      casdm1 = casdm1a + casdm1b
      print(casdm1a)
      print(casdm1b)
      dm1a, dm1b = self.make_full_rdm1s()
      dm1 = dm1a + dm1b
      print(dm1a)
      print(dm1b)
      #print(dm1_beta_mo)
      #--------------------------------------------- nuclear repulsion 
      E_nn   = self.cas._scf.energy_nuc()
      #--------------------------------------------- core parts
      h  = self.cas._scf.get_hcore()
      h_mo = self.ao2mo_transform(self.C_mo,h)
      E_core = self.hcore_energy(h_mo,dm1)
      #--------------------------------------------- J and K parts
      Ja = self.cas._scf.get_j (dm=self.ao2mo_transform(self.C_mo.T,dm1a))
      Jb = self.cas._scf.get_j (dm=self.ao2mo_transform(self.C_mo.T,dm1b))
      Ja_mo = self.ao2mo_transform(self.C_mo,Ja)
      Jb_mo = self.ao2mo_transform(self.C_mo,Jb)
      E_j    = self.Hartree_energy(Ja_mo, Jb_mo, dm1a, dm1b)
      #---------------------------------------------

      print('=======================================================')
      print("   nuclear repulsion energy:     % .8f" % E_nn)
      print("   one-electron energy:          % .8f" % E_core)
      print("   classical Coulomb energy:     % .8f" % E_j)
      #print("   classical Coulomb energy:     % .8f" % E_j1)
      print('=======================================================')

      #casdm1, casdm2 = self.cas.fcisolver.make_rdm1s(self.ci, self.ncas, self.nelecas)
      #casdm1 = self.ao2mo_transform(self.C_mo, casdm1)
      #print(casdm1)
      #print(casdm2)
      #print("\n")
      #print(casdm2)

      return


# ---------------junk----------------
#Te_Vne = np.tensordot (h, dm1)
#E_ja = np.dot(vj_mo, dm1_mo) / 2.0
#vj, vk = cas._scf.get_jk (dm=dm1s)
#vj = vj[0] + vj[1]
#print(np.allclose(dm1, dm1_alpha+dm1_beta))
#dm1s = np.asarray ([dm1_alpha,dm1_beta])

