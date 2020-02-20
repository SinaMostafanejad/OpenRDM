#
# @BEGIN LICENSE
#
# mcpdft by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
from psi4.driver.procrouting import proc

def run_mcpdft(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mcpdft can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mcpdft')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    v2rdm_wfn = kwargs.get('ref_wfn', None)
    if v2rdm_wfn is None:
        raise ValidationError("""Error: %s requires a reference wave function (v2rdm-casscf).""" % name)

    psi4.core.set_variable("V2RDM TOTAL ENERGY",v2rdm_wfn.energy())
   
    if (psi4.core.get_option('MCPDFT', 'MCPDFT_METHOD') == '1DH_MCPDFT'): 
        proc.run_dfmp2('mp2',**kwargs)
    
    func = psi4.core.get_option('MCPDFT','MCPDFT_FUNCTIONAL')
    ref_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    base_wfn = psi4.core.Wavefunction.build(ref_molecule, psi4.core.get_global_option('BASIS'))
    ref_wfn = proc.scf_wavefunction_factory(func, base_wfn, psi4.core.get_option('MCPDFT', 'REFERENCE'))

    # push v2rdm-casscf orbitals onto reference 
    for irrep in range (0,v2rdm_wfn.Ca().nirrep()): 
        ref_wfn.Ca().nph[irrep][:,:] = v2rdm_wfn.Ca().nph[irrep][:,:]
        ref_wfn.Cb().nph[irrep][:,:] = v2rdm_wfn.Cb().nph[irrep][:,:]

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    mcpdft_wfn = psi4.core.plugin('mcpdft.so', ref_wfn)

    return mcpdft_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['mcpdft'] = run_mcpdft


def exampleFN():
    # Your Python code goes here
    pass
