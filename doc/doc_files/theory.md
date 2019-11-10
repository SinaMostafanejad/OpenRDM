Theory     {#theory}
========

[TOC]

@section background Background

Here, we use the conventional notation of multireference (MR) methods
when labeling the molecular orbitals (MOs), \f$\lbrace \psi\rbrace\f$:
the indices \f$i\f$, \f$j\f$, \f$k\f$, and \f$l\f$ denote
inactive (doubly occupied) orbitals; \f$t\f$, \f$u\f$, \f$v\f$, and \f$w\f$ represent
active orbitals; and \f$p\f$, \f$q\f$, \f$r\f$, and \f$s\f$ indicate general orbitals.
[Einstein's summation rule](http://mathworld.wolfram.com/EinsteinSummation.html) is implied in all expressions.

We begin by defining the non-relativistic Born-Oppenheimer (BO) electronic Hamiltonian

\f{equation}{\tag{1}\label{EQ:HAMILTONIAN}
     \hat{H} = h^p_q \hat{E}^p_q + \frac{1}{2} \nu^{pq}_{rs} \hat{e}^{pq}_{rs},
\f}

in terms of one- and two-particle excitation operators which are defined as

\f{gather}{
     \hat{E}^{p}_{q} = \hat{a}^\dagger_{p_\sigma} \hat{a}_{q_\sigma}   \tag{2a}\label{EQ:Epq}  	\\
     \hat{e}^{p r}_{q s} = \hat{E}^p_q \hat{E}^r_s - \delta^q_r \hat{E}^p_s = \hat{a}^\dagger_{p_\sigma} \hat{a}^\dagger_{r_\tau} \hat{a}_{s_\tau} \hat{a}_{q_\sigma} \tag{2b}\label{EQ:epqrs}
\f}

In Eqs. \f$\eqref{EQ:Epq}\f$ and \f$\eqref{EQ:epqrs}\f$, the \f$\hat{a}^\dagger\f$ and \f$\hat{a}\f$ 
represent second-quantized creation and annihilation operators, respectively,
and the Greek labels run over \f$\alpha\f$ and \f$\beta\f$ spins (the sum over which is implied).
The symbol \f$h^p_q = \left<\psi_p|\hat{h}|\psi_q\right>\f$ represents the sum of the electron
kinetic energy and electron-nucleus potential energy integrals, and
\f$\nu^{pq}_{rs} = \left<\psi_p \psi_q|\psi_r\psi_s\right>\f$ is an element of the
two-electron repulsion integral tensor. Because the electronic Hamiltonian
includes up to only pair-wise interactions, the ground-state energy of a
many-electron system can be expressed as an exact linear functional of the
the one-electron reduced-density matrix (1-RDM) and two-electron reduced-density matrix (2-RDM)

\f{equation}{\tag{3}\label{EQ:Eel}
E = {}^1D^p_q h^p_q + \frac{1}{2} {}^2D^{pq}_{rs} \nu^{pq}_{rs}.
\f}

Here, the 1-RDM and the 2-RDM are represented in their spin-free forms, 
with elements defined as

\f{equation}{\tag{4a}
{}^1D^p_q = {}^1D^{p_\sigma}_{q_\sigma} = \left<\Psi|\hat{a}^\dagger_{p_\sigma} \hat{a}_{q_\sigma}|\Psi\right>	\label{EQ:1RDM}
\f}

and

\f{equation}{\tag{4b}
{}^2D^{pq}_{rs} = {}^2D^{p_\sigma q_\tau}_{r_\sigma s_\tau} = \left<\Psi|\hat{a}^\dagger_{p_\sigma} \hat{a}^\dagger_{q_\tau} \hat{a}_{s_\tau} \hat{a}_{r_\sigma}|\Psi\right> \label{EQ:2RDM},
\f}

respectively. Again, the summation over the spin labels 
in Eqs. \f$\eqref{EQ:1RDM}\f$ and \f$\eqref{EQ:2RDM}\f$ is implied.

@section mcpdft Multiconfiguration Pair-Density Functional Theory


The MCDPFT energy expression can be written as

\f{equation}{\tag{5}\label{EQ:EMCPDFT}
E_{\text{MCPDFT}} = 2h^i_i + h^t_u {}^1D^t_u + E_\text{H} + E_\text{xc}\left[\rho,\Pi,|\nabla\rho|,|\nabla\Pi|\right],
\f}

where the Hartree energy, \f$E_\text{H}\f$, is

\f{equation}{\tag{6}\label{EQ:EHARTREE}
E_\text{H} = 2 \nu^{ij}_{ij} + 2\nu^{ti}_{ui} {}^1D^t_u + \frac{1}{2} \nu^{tv}_{uw} {}^1D^{t}_{u} {}^1D^{v}_{w}
\f}

The total electronic density and its gradient are defined by the 1-RDM as

\f{equation}{\label{EQ:RHO}\tag{7}
\rho(\mathbf{r}) = {}^1D^p_q\ \psi^*_p(\mathbf{r}) \psi_q(\mathbf{r}),
\f}

and

\f{equation}{\label{EQ:DRHO}\tag{8}
\nabla\rho(\mathbf{r}) = {}^1D^p_q \left[ \nabla\psi^*_p(\mathbf{r}) \psi_q(\mathbf{r}) + \psi^*_p(\mathbf{r}) \nabla\psi_q(\mathbf{r}) \right],
\f}

respectively. The on-top pair density (OTPD) and its gradient can similarly be defined in
terms of the 2-RDM as

\f{equation}{
\label{EQ:PI}\tag{9}
\Pi(\mathbf{r})  = {}^2D^{pq}_{rs}\ \psi^*_p(\mathbf{r}) \psi^*_q(\mathbf{r}) \psi_r(\mathbf{r}) \psi_s(\mathbf{r}),
\f}

and

\f{eqnarray}{
\label{EQ:DPI}\tag{10}
\nabla\Pi(\mathbf{r}) = {}^2D^{pq}_{rs} &[& \nabla\psi^*_p(\mathbf{r}) \psi^*_q(\mathbf{r}) \psi_r(\mathbf{r}) \psi_s(\mathbf{r}) \nonumber \\\
                                 &+& \psi^*_p(\mathbf{r}) \nabla\psi^*_q(\mathbf{r}) \psi_r(\mathbf{r}) \psi_s(\mathbf{r}) \nonumber \\\
                                 &+& \psi^*_p(\mathbf{r}) \psi^*_q(\mathbf{r}) \nabla\psi_r(\mathbf{r}) \psi_s(\mathbf{r}) \nonumber \\\
                                 &+& \psi^*_p(\mathbf{r}) \psi^*_q(\mathbf{r}) \psi_r(\mathbf{r}) \nabla\psi_s(\mathbf{r}) ~],
\f}

respectively. Here, the 1- and 2-RDMs are obtained from an MR computation.
