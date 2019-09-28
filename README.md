
<p align="center">
<!--img src="logo.png" style='height: 30%; width: 50%; object-fit: contain'/--> 
<br>
</p>

<table align="center">
  <tr>
     <th align="left">CI/CD</th>
     <th align="left">
        <a href="https://travis-ci.com/SinaMostafanejad/OpenRDM"><img alt="Travis CI" src="https://travis-ci.com/SinaMostafanejad/OpenRDM.svg?token=aVpZaqKz4Vv5czxgJ8WE&branch=master"></a>
        <a href="https://ci.appveyor.com/project/SinaMostafanejad/openrdm"><img alt="AppVeyor" src="https://ci.appveyor.com/api/projects/status/67t0souy2fhoc7l5?svg=true"></a>
     </th>
  </tr>
  <tr>
    <th align="left">Code Coverage and Quality</th>
    <th align="left">
      <a href="https://codecov.io/gh/SinaMostafanejad/OpenRDM">
<img alt="CodeCoverage" src="https://codecov.io/gh/SinaMostafanejad/OpenRDM/branch/master/graph/badge.svg" />
      </a>
      <a href="https://lgtm.com/projects/g/SinaMostafanejad/OpenRDM/context:cpp"><img alt="Language grade: C/C++"       src="https://img.shields.io/lgtm/grade/cpp/g/SinaMostafanejad/OpenRDM.svg?logo=lgtm&logoWidth=18"/></a> 
<a href="https://www.codefactor.io/repository/github/sinamostafanejad/openrdm"><img src="https://www.codefactor.io/repository/github/sinamostafanejad/openrdm/badge" alt="CodeFactor" /></a>
     </th>
  </tr>
  <tr>
    <th align="left">Foundation</th>
    <th align="left">
      <a href="https://opensource.org/licenses/BSD-3-Clause"><img alt="GitHub license" src="https://img.shields.io/badge/license-BSD--3-blueviolet"></a>
      <a href="https://www.linuxfoundation.org/"><img alt="GitHub license" src="https://img.shields.io/badge/Platforms-Linux-blue"></a>
    </th>
  </tr>
  <tr>
    <th align="left">Miscellaneous</th>
    <th align="left">
      <a href="http://hits.dwyl.io/SinaMostafanejad/OpenRDM"><img alt="Hit counts" src="http://hits.dwyl.io/SinaMostafanejad/OpenRDM.svg"></a>
    </th>
  </tr>
</table>


# OpenRDM

An ab initio library for strongly-correlated many-body systems based on multiconfiguration pair-density functional theory.

## OVERVIEW

<b>OpenRDM</b> is a source-independent version of the original <a href="https://github.com/edeprince3/RDMinoles">RDM-INOLES</a> plugin to the Psi4 quantum chemistry program package. <b>OpenRDM</b> uses multiconfigurational pair-density functional theory (MCPDFT) to provide an accurate and efficient description of static and dynamical correlation effects.

Please refer to the <b>OpenRDM</b>' <a href="https://sinamostafanejad.github.io/libRDMInoles/index.html">documentation</a> for further details about MCPDFT and its implementation.
<!-- Both translated and fully-translated versions of Slater and Vosko-Wilk-Nusair random-phase approximation expression III (SVWN3), Perdew-Burke-Ernzerhof (PBE), revised PBE (revPBE), Becke88 exchange and one-parameter correlation functional (BOP) and Becke and Lee-Yang-Parr (BLYP) on-top pair-density exchange-correlation functionals are available at the moment. In addition, the global-, double- and range-separated hybrid multi-configurational OTPDs such as wPBE and LRC-wPBE have also been implemented. However, this part of the project also is under the ongoing developement.

In summary, RDM-INOLES:

* can provide an interface with any (multiconfigurational) method that is able to provide 1-electron and 2-electron RDMs.
* hosts the variational 2-RDM driven complete active-space self-consistent field (v2RDM-CASSCF) as the reference method [2] by default
* can generate a .wfn file for further analysis of the wavefunction based on the quantum theory of atoms in molecules (QTAIMs)
* uses the reference total density and on-top pair-density (OTPD) functions as the input to build the so-called OTPD exchange-correlation (XC) functionals [1]
* features a double-hybrid MCPDFT method that is based on the linearly-scaled one-parameter double-hybrid (LS1DH) of Toulouse et al. described in Ref [3]
* will include E. Valeev's universal perturbative explicitly correlated basis-set incompleteness correction [4]
* will provide and support both scaled and unscaled densities in MCPDFT
-->

## INSTALLATION

We adopt CMake for the installation procedure and package management.

At the moment, this prototype can be processed by running the configure script to obtain an MCPDFT energy correction on the H2 molecule at its equilibrium bond length adopting minimal basis and SVWN3 exchange-correlation functional.

Please feel free to modify the configure script and/or CMakeLists.txt to address your needs.

## REFERENCES

[1] M. Mostafanejad and A. E. DePrince III, J. Chem. Theory Comput. 15, 290-302 (2019). "Combining Pair-Density Functional Theory and Variational Two-Electron Reduced-Density Matrix Methods"

[2] J. Fosso-Tande, T.-S. Nguyen, G. Gidofalvi, and A. E. DePrince III, J. Chem. Theory Comput., 12, 2260-2271 (2016). "Large-scale variational two-electron reduced-density-matrix-driven complete active space self-consistent field methods."

[3] J. Toulouse, K. Sharkas, E. Bremond and C. Adamo J. Chem. Phys. 135, 101102 (2011). "Rationale for a new class of double-hybrid approximations in density-functional theory"

[4] M. Torheyden and E. F. Valeev J. Chem. Phys. 131, 171103 (2009). "Universal perturbative explicitly correlatedbasis set incompleteness correction"
