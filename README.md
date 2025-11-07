<!-- Title -->
<h1 align="center">
  FLOW
</h1>

<!-- Information badges -->
<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg" />
  </a>
  <a href="https://fortran-lang.org">
    <img alt="Fortran" src="https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat" />
  </a>
</p>

<p align="center">
    <img width="300" src="https://github.com/bradyelster/FLOW/blob/046650ccd7100a337e5b64f068555b5b476c0e1e/docs/transonic_rho.png">
</p>

FLOW was developed as a collaboration between the University of Rochester and the Princeton Plasma Physics Laboratory. FLOW is written primarily in Fortran and is used to study the equilibrium properties of toroidal devices, such as tokamaks, in conditions relevant to present day experiments. 

The most unique feature of the code is the ability to study flow-dependent equilibria. Present day experiments often show the presence of high macroscopic flow. Large toroidal flows, of the order of a significant fraction of the sound speed, are routinely measured in tokamak plasmas (e.g., DIII-D, Alcator C-Mod), even in the absence of an external source of momentum. Measured toroidal rotation is particularly large in low-aspect ratio machines (NSTX, MAST). Poloidal rotation in experiments is somewhat lower and mostly observed near the plasma edge. Edge rotation can be of the order of the poloidal sound speed (~10s km/s): even such a small rotation can create large qualitative modifications in the equilibrium profiles. [R. Betti, J. P. Freidberg, PoP 7, 2439 (2000)](https://pubs.aip.org/aip/pop/article-abstract/7/6/2439/103410/Radial-discontinuities-in-tokamak?redirectedFrom=fulltext)

Despite the obvious relevance of the study of equilibrium in the presence of strong flow, little work has been done on the subject, in particular regarding poloidal rotation, which constitutes a more difficult problem with respect to purely toroidal rotation. To the best of the authors' knowledge, before the code development no numerical work at all had been done on equilibria with poloidal flow exceeding the poloidal sound speed.

The code main targets are:
1. The study of MHD and kinetic tokamak equilibria with purely toroidal flow.
2. The study of MHD equilibria in the presence of arbitrary poloidal flow.

For a more detailed description of the code, see [L. Guazzotto, R. Betti, J. Manickam and S. Kaye, PoP 11, 604 (2004)](https://pubs.aip.org/aip/pop/article-abstract/11/2/604/260879/Numerical-study-of-tokamak-equilibria-with?redirectedFrom=fulltext)

## Links

* Project homepage: <https://www.auburn.edu/cosam/faculty/physics/guazzotto/research/FLOW_main.html>
* Author homepage: <https://www.auburn.edu/cosam/faculty/physics/guazzotto/research/>
* Code repository: <https://github.com/bradyelster/FLOW>
* Documentation: <http://FLOW-framework.readthedocs.org>
