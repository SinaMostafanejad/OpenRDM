import numpy as np
import pyscf
import h5py
from scipy import sparse

from pyscf import gto, dft, scf, fci
from pyscf.dft import numint
#from pyscf.tools import cubegen

mol = pyscf.M(
    atom = 'H 0 0 0; H 0 0 0.75',  # in Angstrom
    basis = 'sto3g',
    symmetry = True,
)

#myhf = mol.HF()
#myhf.kernel()
myhf = dft.RKS(mol)
myhf.grids.atom_grid = (75, 302)

#myhf.grids.radi_method = dft.treutler
myhf.grids.radi_method = dft.mura_knowles

myhf.grids.becke_scheme = dft.original_becke

myhf.grids.prune = dft.nwchem_prune

#myhf.grids.level = 3

myhf.kernel()

# Orbital energies, Mulliken population etc.
#myhf.analyze()

nao = mol.nao_nr()

cmat_mo = myhf.mo_coeff
print("AO2MO transformation matrix C =")
print(cmat_mo)
cmat_mo_occ = myhf.mo_occ

smat = myhf.get_ovlp()
print("overlap matrix S =")
print(smat)

#dmat = myhf.make_rdm1()
#dtemp = np.matmul(cmat_mo.transpose(),dmat)
#dmat_mo = np.matmul(dtemp,cmat_mo)
dmat_mo = np.zeros((nao,nao))
dmat_mo[0][0] = 1.0
print(dmat_mo)
dm = myhf.make_rdm1()

print(dm)

#eigval, eigvec = myhf.eig(dmat_mo,smat)
#print(eigval)

hcore_ao = myhf.get_hcore()
htemp = np.matmul(cmat_mo.transpose(),hcore_ao)
hcore_mo = np.matmul(htemp,cmat_mo)
print(hcore_mo)

core_en = 2.0 * np.vdot(dmat_mo,hcore_mo)
print("core energy = %15.10lf\n" % core_en)

jmat_ao = myhf.get_j()
jtemp = np.matmul(cmat_mo.transpose(),jmat_ao)
jmat_mo = np.matmul(jtemp,cmat_mo)
jmat_mo *= 0.5
print(jmat_mo)

coulomb_en = 2.0 * np.vdot(dmat_mo,jmat_mo)
print("Classical Coulomb energy = %15.10lf\n" % coulomb_en)

coords = myhf.grids.coords
weights = myhf.grids.weights
ao_value = numint.eval_ao(mol, coords, deriv=1)
f = h5py.File("dataset.h5",'w')
f["/Grids"] = ao_value
rho0 = f["/Grids"][0,:]
print(rho0)
print(ao_value.shape)
#print(coords.shape)
# The first row of rho is electron density, the rest three rows are electron
# density gradients which are needed for GGA functional
rho = numint.eval_rho(mol, ao_value, dmat_mo, xctype='GGA')

