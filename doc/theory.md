#### Theory          {#theory}

Here, we use the conventional notation of multireference (MR) methods
when labeling the molecular orbitals (MOs), \f$\lbrace \psi\rbrace\f$:
the indices \f$i\f$, \f$j\f$, \f$k\f$, and \f$l\f$ denote
inactive (doubly occupied) orbitals; \f$t\f$, \f$u\f$, \f$v\f$, and \f$w\f$ represent
active orbitals; and \f$p\f$, \f$q\f$, \f$r\f$, and \f$s\f$ indicate general orbitals. A
summation over repeated indices is implied in all expressions.

We begin by defining the non-relativistic Born-Oppenheimer (BO) electronic Hamiltonian

\begin{equation}
     \hat{H} = h^p_q \hat{E}^p_q + \frac{1}{2} \nu^{pq}_{rs} \hat{e}^{pq}_{rs}
\end{equation}

in which, the one- and two-particle excitation operators can be expressed as

\f{gather}{
     \hat{E}^{p}_{q} = \hat{a}^\dagger_{p_\sigma} \hat{a}_{q_\sigma}   \label{EQ:Epq}  	\\
     \hat{e}^{p r}_{q s} = \hat{E}^p_q \hat{E}^r_s - \delta^q_r \hat{E}^p_s = \hat{a}^\dagger_{p_\sigma} \hat{a}^\dagger_{r_\tau} \hat{a}_{s_\tau} \hat{a}_{q_\sigma} \label{EQ:epqrs}
\f}

where \f$\hat{a}^\dagger\f$ and \f$\hat{a}\f$ represent second-quantized creation and
annihilation operators, respectively, and the Greek labels run over \f$\alpha\f$
and \f$\beta\f$ spins (the sum over which is implied). The symbol 


The MCDPFT energy expression can be written as

\begin{equation}\label{EQ:EMCPDFT}
E_{\text{MCPDFT}} = 2h^i_i + h^t_u {}^1D^t_u + E_\text{H} + E_\text{xc}\left[\rho,\Pi,|\nabla\rho|,|\nabla\Pi|\right],
\end{equation}

where Eq. \f$\ref{EQ:epqrs}\f$
