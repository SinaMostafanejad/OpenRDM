import h5py
import numpy as np
import pyscf

from pyscf import gto, scf, fci
#from pyscf.tools import cubegen

mol = pyscf.M(
    atom = 'H 0 0 0; H 0 0 0.75',  # in Angstrom
    basis = 'sto3g',
    symmetry = True,
)

myhf = mol.HF()
myhf.kernel()

# Orbital energies, Mulliken population etc.
#myhf.analyze()

nao = mol.nao_nr()

cmat_mo = myhf.mo_coeff
print(cmat_mo)
cmat_mo_occ = myhf.mo_occ

overlap = myhf.get_ovlp()

#dmat = myhf.make_rdm1()
#dtemp = np.matmul(cmat_mo.transpose(),dmat)
#dmat_mo = np.matmul(dtemp,cmat_mo)
dmat_mo = np.zeros((nao,nao))
dmat_mo[0][0] = 1.0
print(dmat_mo)

eigval, eigvec = myhf.eig(dmat_mo,overlap)
print(eigval)

jmat_ao = myhf.get_j()
jtemp = np.matmul(cmat_mo.transpose(),jmat_ao)
jmat_mo = np.matmul(jtemp,cmat_mo)
jmat_mo *= 0.5
print(jmat_mo)

ca = 2.0 * np.vdot(dmat_mo,jmat_mo)
print("Coulomb energy = %15.10lf\n" % ca)


