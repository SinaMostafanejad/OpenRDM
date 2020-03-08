import numpy as np
import scipy as sp
from functools import reduce
import time
import h5py
import pyscf
#from pyscf import dft, ao2mo, fci, mcscf
#from pyscf.lib import logger, temporary_env
#from pyscf.mcscf import mc_ao2mo
#from pyscf.mcscf.addons import StateAverageMCSCFSolver
from pyscf import gto, scf, mcscf

bas_lst = ['sto-3g','ccpvdz']

mol = gto.M(
    atom = 'H 0 0 0; H 0 0 0.75',
    basis = bas_lst[0])

hf = scf.RHF(mol).run()

# CAS(2e,2o)
nele = 2
norb = 2
cas = mcscf.DFCASSCF(hf, norb, nele)
cas.natorb = True
cas.kernel()

#===================================================
C_mo = cas.mo_coeff
ci = cas.ci
nao = C_mo.shape[0]
nmo = C_mo.shape[1]
nfrz = cas.frozen
ncore = cas.ncore
ncas = cas.ncas
nelecas = cas.nelecas
#---------------------------------------------------
#casdm1, casdm2 = cas.fcisolver.make_rdm12(ci, ncas, nelecas)
#print(casdm1)
#print("\n")
#print(casdm2)
#---------------------------------------------------
print("\n")
#print("Frozen orbitals:     %d" % nfrz)
print("Number of AOs:     %s" % nao)
print("Number of MOs:     %s" % nmo)
print("Inactive (core) space size:     %s" % ncore)
print("Active space size:     %s" % ncas)
print("Active(alpha,beta) electrons:     (%s,%s)" % nelecas)
print("\n")
#===================================================
# 1pdm in AO representation
dm1 = cas.make_rdm1()
# alpha and beta 1-pdm in AO representation
dm1_alpha, dm1_beta = cas.make_rdm1s()
print(np.allclose(dm1, dm1_alpha+dm1_beta))
dm1s = np.asarray ([dm1_alpha,dm1_beta])
#---------------------------------------------------
#dm1_temp = np.matmul(C_mo.transpose(),dm1)
#dm1_mo = np.matmul(dm1_temp,C_mo)
#print(dm1_mo)
#print("\n")
#---------------------------------------------------
dm1_mo = reduce(np.dot,(C_mo.T,dm1,C_mo))
#print(dm1_mo)
#print("\n")
#===================================================
#vj, vk = cas._scf.get_jk (dm=dm1s)
#vj = vj[0] + vj[1]
#print(vj)
Ja = cas._scf.get_j (dm=dm1_alpha)
Jb = cas._scf.get_j (dm=dm1_beta)
Ja_mo = reduce(np.dot,(C_mo.T,Ja,C_mo))
Jb_mo = reduce(np.dot,(C_mo.T,Jb,C_mo))
print(Ja)
print(Jb)
Ej_aa = reduce(np.tensordot,(Ja,dm1_alpha))
Ej_ab = reduce(np.tensordot,(Ja,dm1_beta))
Ej_ba = reduce(np.tensordot,(Jb,dm1_alpha))
Ej_bb = reduce(np.tensordot,(Jb,dm1_beta))
E_j = 0.5 * ( Ej_aa + Ej_ab + Ej_ba + Ej_bb )
print("classical Coulomb energy:     %s" % E_j)
#===================================================
#Te_Vne = np.tensordot (h, dm1)
Vnn = cas._scf.energy_nuc()
h = cas._scf.get_hcore()
print()
#E_ja = np.dot(vj_mo, dm1_mo) / 2.0

