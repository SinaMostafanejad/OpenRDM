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

cmat = myhf.mo_coeff
c_occ = myhf.mo_occ

overlap = myhf.get_ovlp()

dmat = myhf.make_rdm1()
dtemp = np.matmul(cmat.transpose(),dmat)
dmat_mo = np.matmul(dtemp,cmat)
print(dmat_mo)

eigval, eigvec = myhf.eig(dmat,overlap)
print(eigval)

jmat_ao = myhf.get_j()
jtemp = np.matmul(cmat.transpose(),jmat_ao)
jmat_mo = np.matmul(jtemp,cmat)
print(0.5 * jmat_mo)
