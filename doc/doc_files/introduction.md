# Introduction         {#introduction}

The <b>OpenRDM</b> library is an open-source platform for reduced-density matrix-based analysis and computation.
It adopts the multiconfiguration pair-density functional theory (MCPDFT) to provide an accurate and efficient
description of the electronic structure of systems with strong multireference (MR) character.

The accurate and computationally affordable description of the electronic
structure of many-body systems remains a major challenge within the
quantum chemistry and molecular physics
communities. Specifically, the realization of general approaches that account for both
dynamical and nondynamical correlation effects in multiconfigurational or strongly-correlated
systems is particularly difficult. The main issue is that many approaches designed to deal with MR
problems are not particularly efficient for capturing dynamical correlation
effects. A similar statement can be made regarding the ability of methods
designed to model dynamical correlation to capture MR effects.

One can broadly classify approaches to the electron correlation problem as
either falling within wave function theory (WFT), in which the many-electron wave function
is the central quantity, or density-based theories,
which include both density functional theory (DFT) and
RDM-based approaches. In principle, WFT is preferable, as it allows for systematic
improvement in the calculated energies and properties of the
system. In practice, however, the computational complexity of post-HF wave-function-based methods,
specifically MR approaches, limits their application to small systems.
The wide-ranging success of DFT, on the other hand, stems from the its ability to provide
a reasonable description of electron correlation at significantly lower costs.  Nonetheless, DFT
often fails for MR systems, and it does not offer a systematic approach for improving its accuracy.
