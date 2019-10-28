Introduction         {#introduction}
==========

The <b>OpenRDM</b> library is an open-source platform for reduced-density matrix-based analysis and computation.
It adopts the multiconfiguration pair-density functional theory (MCPDFT) to provide an accurate and efficient
description of the electronic structure of systems with strong multireference (MR) character. The MCPDFT addresses
the double counting in the (Coulomb) correlation and the symmetry dilemma through a change of variables in the 
Kohn-Sham DFT (KSDFT) exchange-correlation (XC) functionals from spin-densities to total and on-top pair densities.
However, the MCPDFT leaves the problem of the computational cost barrier of the MR methods open.

In order to address the poor scaling of the conventional configuration-interaction (CI)-based methods
within the framework of MCPDFT, we combine the variational 2-electron reduced-density matrix (v2RDM)-driven
complete active-space self-consistent field (CASSCF) method with pair-density functional theory.
The v2RDM-CASSCF represents the electronic structure of the active space with the 2-RDM as opposed to 
the CI wave fucntion leading to the realization of the polynomial cost with respect to the active space size.
As such, the v2RDM-CASSCF method facilitates the v2RDM-based MCPDFT calculations with active spaces as large
as 64 electrons in 64 orbitals. For our recent studies on strongly correlated systems with large active spaces
see [here](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00988)
and [here](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00768).

In the [Theory](@ref theory) section, we delve into the details of the MCPDFT. In order to see <b>OpenRDM</b>
in action, take a glance at the [How to Use](@ref howtouse) section.
