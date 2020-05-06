import numpy as np
import scipy as sp
from functools import reduce
import time
import h5py
import pyscf
from pyscf import gto, scf, mcscf
from pyscf.dmrgscf import dmrgci

#---------------------------------------------------
import libpyscf
#---------------------------------------------------

bas_lst    = ['sto-3g','ccpvdz']
method_lst = ['MCSCF', 'DMRG']

mol = gto.M(
    atom = 'H 0 0 0; H 0 0 0.75',
    basis = bas_lst[1])

hf = scf.RHF(mol).run()

# CAS(2e,2o)
nele = 2
norb = 4
#===================================================
#                   MCSCF + PDFT
#=================================================== MCSCF
cas = mcscf.DFCASSCF(hf, norb, nele, auxbasis='def2-svp-jkfit')     # w/ density-fitting
#cas = mcscf.CASSCF(hf, norb, nele)                                 # w/o density fitting
#cas.natorb = True
cas.kernel()
#--------------------------------------------------- PDFT
mcpdft = libpyscf.MCPDFT(mol, cas, ref_method=method_lst[0])
mcpdft.write_hdf5 = True
mcpdft.kernel()
#===================================================

#===================================================
#                   DMRGSCF + PDFT
#=================================================== DMRGSCF
cas = dmrgci.DMRGSCF(hf, norb, nele, maxM=2000, tol=1.e-8)
cas.kernel()
#--------------------------------------------------- PDFT
mcpdft = libpyscf.MCPDFT(mol, cas, ref_method=method_lst[1])
mcpdft.write_hdf5 = True
mcpdft.kernel()
#===================================================
