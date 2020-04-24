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
#---------------------------------------------------
import libpyscf
#---------------------------------------------------

bas_lst = ['sto-3g','ccpvdz']

mol = gto.M(
    atom = 'H 0 0 0; H 0 0 0.75',
    basis = bas_lst[0])

hf = scf.RHF(mol).run()

# CAS(2e,2o)
nele = 2
norb = 2
cas = mcscf.DFCASSCF(hf, norb, nele, auxbasis='def2-svp-jkfit')
#cas = mcscf.CASSCF(hf, norb, nele, auxbasis='def2-svp-jkfit')
cas.natorb = True
cas.kernel()
#===================================================
mcpdft = libpyscf.MCPDFT(cas)
mcpdft.write_hdf5 = True
mcpdft.kernel(mol)
#===================================================
