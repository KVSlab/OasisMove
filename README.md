# OasisMove - Moving Domain CFD Solver

_________________
[![GPL-3.0](https://img.shields.io/github/license/kvslab/oasismove)](LICENSE)
[![codecov](https://codecov.io/gh/KVSlab/OasisMove/branch/main/graph/badge.svg?token=ETTVXYFHJ2)](https://codecov.io/gh/KVSlab/OasisMove)
[![CI](https://github.com/kvslab/oasismove/actions/workflows/check_and_test_package.yml/badge.svg)](https://github.com/kvslab/oasismove/actions/workflows/check_and_test_package.yml)
[![GitHub pages](https://github.com/kvslab/oasismove/actions/workflows/deploy_pages.yml/badge.svg)](https://github.com/kvslab/oasismove/actions/workflows/deploy_pages.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7868226.svg)](https://doi.org/10.5281/zenodo.7868226)
_________________

<p align="center">
    <img src="docs/figures/moving_atrium.gif" width="640" height="315" alt="Left atrium flow"/>
</p>
<p align="center">
    Fluid velocity (left) and mesh deformation (right) of a moving patient-specific left atrium model, simulated over two cardiac cycles lasting for $T=2$ s.
    The model is publicly available from <a href="https://doi.org/10.5281/zenodo.5801337">this</a> dataset. 
</p>

Description
-----------
OasisMove is a high-level/high-performance open-source Navier-Stokes solver for fluid flow in rigid and moving domains
written in Python/[FEniCS](https://fenicsproject.org/), and is an extension of the computational fluid dynamics (CFD)
solver [Oasis](https://github.com/mikaem/Oasis). In OasisMove, the Navier-Stokes equations are expressed in the arbitrary
Lagrangian-Eulerian formulation, which is suitable for handling moving domains. This moving domain solver has undergone
rigorous verification and validation, and results have shown that OasisMove follows theoretical convergence rates, being
second order accurate in time, and second and third order accurate in space with P1/P1 and P2/P1 finite elements.
OasisMove has been developed with cardiovascular flows in mind, but is applicable to several flow problems within CFD.

<p align="center">
    <img src=docs/figures/verification_u_p.png width="630 height="470" alt="Convergence rate analysis"/>
</p>
<p align="center">
    Spatial convergence study of OasisMove performed by varying the characteristic edge length Δx. On the left, the L2 error for the
    velocity, and on the right the L2 error for the pressure, both following theoretical convergence rates. 
    The solid lines represent the simulation results, and the dashed lines display the theoretical convergence rates. 
    A similar study was performed to address temporal convergence, resulting in second order convergence (not shown here). 
</p>


Installation
------------
OasisMove and its dependencies can be installed using either `conda`, or by building and running a `Docker` container,
and `pip`. For detailed installation notes see
the [installation guidelines](https://kvslab.github.io/OasisMove/installation.html).

Documentation
-------------
OasisMove's documentation is hosted [here](https://kvslab.github.io/OasisMove). This includes
multiple [tutorials](https://kvslab.github.io/OasisMove/tutorials.html), meant to guide the user through the basic steps
of performing a computational fluid dynamic simulation and creating problem files.

For futher details on vanilla Oasis, please refer to its [wiki](https://github.com/mikaem/oasis/wiki) or
the [user manual](https://github.com/mikaem/Oasis/tree/master/doc/usermanual.pdf)

If you wish to use OasisMove for journal publications, please cite the
following [paper](https://onlinelibrary.wiley.com/doi/10.1002/cnm.3703).

Licence
-------
OasisMove is licensed under the GNU GPL, version 3 or (at your option) any later version.

OasisMove is Copyright (2018-2023) by the authors.

Authors
-------
OasisMove has been developed by

* [Henrik A. Kjeldsberg](https://github.com/HKjeldsberg)

Issues
------
Please report bugs and other issues through the issue tracker at:

https://github.com/KVSlab/OasisMove/issues
