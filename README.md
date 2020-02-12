
<p align="center">
<img src="OpenRDM.png" style='height: 1%; width: 5%; object-fit: contain'/> 
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
    <th align="left">Technical Support</th>
    <th align="left">
      <a href="https://join.slack.com/t/openrdm/shared_invite/enQtOTM2MDg2MzUxNjIyLWNlMzFlOWFhYTVlZGQ3ZGYxNWY3NTk4ZjRhYzM3MTU5MWZhN2VhY2Y5NzBiNjVjYzU1YWJkZDc2ODdhYTM4Yjg"><img alt="Slack Link" src="https://img.shields.io/badge/Chat on-Slack-blue?style=flat&logo=slack"></a>
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
     <!-- <a href="http://hits.dwyl.io/SinaMostafanejad/OpenRDM"><img alt="Hit counts" src="http://hits.dwyl.io/SinaMostafanejad/OpenRDM.svg"></a> -->
Visitors: <a href="https://github.com/SinaMostafanejad/OpenRDM" ><img src="https://www.webfreecounter.com/hit.php?id=guxcaqf&nd=9&style=20" border="0" alt="free counter"></a>
    </th>
  </tr>
</table>


# OpenRDM

An open-source library for reduced-density matrix-based analysis and computation.

## OVERVIEW

<b>OpenRDM</b> is a standalone quantum chemistry software that adopts multiconfigurational pair-density functional theory (MCPDFT) to provide an accurate and efficient description of static and dynamical correlation effects in strongly-correlated systems.

Its plugin to the Psi4 quantum chemistry program package can be found <a href="https://github.com/edeprince3/RDMinoles">here</a>.

Please refer to the <a href="https://sinamostafanejad.github.io/OpenRDM/index.html">documentation</a> for further details about <b>OpenRDM</b>.
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

The installation procedure is detailed in the [documentation](https://sinamostafanejad.github.io/OpenRDM/howtouse.html).

For user's convenience, we have placed a configure script in project's root directory; Feel free to change its contents to address your needs.

## KNOWN ISSUES

If you have any problems working with <b>OpenRDM</b> software, please take a glance at the [Known Issues](https://sinamostafanejad.github.io/OpenRDM/issues.html) page of the [documentation](https://sinamostafanejad.github.io/OpenRDM/index.html).

You can report bugs or request new features by opening issues and submitting pull requests. Feel free to contact the author for any feedbacks or suggestions.

## HOW TO CITE US

We appreciate your time for using <b>OpenRDM</b> and do our best to keep the program up-to-date. Your support through citing the following manuscripts helps us to add more exciting features to <b>OpenRDM</b> and improve its performance in its future releases.

[1] [M. Mostafanejad and A. E. DePrince III, J. Chem. Theory Comput. 15, 290-302 (2019). "Combining Pair-Density Functional Theory and Variational Two-Electron Reduced-Density Matrix Methods"](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00988)

[2] [M. Mostafanejad, M. D. Liebenthal, and A. E. DePrince III arXiv:1911.11162 (2019). "Global hybrid multiconfiguration pair-density functional theory"](https://arxiv.org/abs/1911.11162v1)

[3] [J. W. Mullinax, L. N. Koulias, E. Maradzike, M. Mostafanejad, E. Epifanovsky, G. Gidofalvi, and A. E. DePrince III, J. Chem. Theory Comput. 15, 6164-6178 (2019). "Heterogeneous CPU + GPU Algorithm for Variational Two-Electron Reduced-Density Matrix-Driven Complete Active-Space Self-Consistent Field Theory"](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00768)

[4] [J. Fosso-Tande, T.-S. Nguyen, G. Gidofalvi, and A. E. DePrince III, J. Chem. Theory Comput., 12, 2260-2271 (2016). "Large-scale variational two-electron reduced-density-matrix-driven complete active space self-consistent field methods"](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00190)
