---
title: 'OasisMove: A Verified and Validated Computational Fluid Dynamics Solver for Moving Domains'
tags:
- Python
- computational fluid dynamics
- FEniCS
- movind domain
- simulation

authors:
- name: Henrik A. Kjeldsberg
  orcid: 0000-0002-7764-4248
  affiliation: 1
- name: Joakim Sundnes
  orcid: 0000-0002-1890-7722
  affiliation: 1
- name: Kristian Valen-Sendstad
  orcid: 0000-0002-2907-0171
  affiliation: 1

affiliations:
- name: Department of Computational Physiology, Simula Research Laboratory
  index: 1

date: 31 May 2023
bibliography: paper.bib
---

# Summary

Computational fluid dynamics (CFD) has emerged as a valuable alternative to costly and time-consuming laboratory
experiments in various scientific and engineering disciplines. The intricate interaction between fluid dynamics and
physical processes, relevant in biomedical and cardiovascular applications, presents an ongoing challenge in CFD for
achieving accurate and efficient representations of fluid flow in moving domains. By incorporating the Arbitrary
Lagrangian-Eulerian (ALE) formulation of the Navier-Stokes equations, `OasisMove` addresses the added complexity of
moving domains, allowing for seamless fluid flow simulations in both rigid and moving domains. The software is
implemented in the finite element method framework FEniCS, and provides an easy learning curve and wide applicability.
Rigorous verification and validation procedures ensure the dependability and precision of `OasisMove`, displaying
excellent agreement with theoretical convergence rates, numerical benchmarks, and experimental results.

As a general open-source solver with entry level high-performance computing capabilities, `OasisMove` is well-suited for
addressing diverse fluid challenges. Its robust performance and proven capabilities in handling transitional
cardiovascular flow problems make it a valuable asset for researchers and engineers seeking reliable and adaptable fluid
dynamics solutions.

# Statement of Need

Over the past decades, CFD simulations using patient-specific medical images have shown promising correlations between
flow phenotype and disease initiation, progression and
outcome [@taylor2010image; @taylor2013computational; @fisher; @Lee2005a; @steinman2003image]. However, the variety of
numerical solution strategies used in the literature can lead to inconsistencies in results, and make direct comparisons
between studies difficult. While commercial solvers such as Ansys Fluent [@matsson2021introduction] or
COMSOL [@pryor2009multiphysics] offer user-friendly interfaces and are relatively accessible for non-experts, their
default settings can sacrifice accuracy for speed [@valen2014mind], and users have limited control over the underlying
numerical methods and algorithms employed. On the other hand, open-source solvers like OpenFOAM [@jasak2007openfoam] or
LifeV [@bertagna2017lifev] offer increased flexibility but necessitate a more profound understanding of numerical
methods and present a steeper learning curve for users. Additionally, there is notable variability in solution
strategy [@valen2018real], partly due to factors such as unavailable source codes or open-source solvers implemented in
extensive and overwhelming low-level libraries. The observed variability highlights the need for increased transparency.

The aforementioned challenges are further exacerbated by the added feature of moving domains, which is important for
many cardiovascular flow problems. The common assumption of rigid domains has been successful for applications with
minimal deformations, such as smaller blood vessels [@valen2014high], but it is generally not applicable for CFD
applications in domains such as the heart and major blood
vessels  [@jia2019image; @sanatkhani2021subject; @pasta2022inversion; @sanatkhani10subject]. In these scenarios,
assuming rigid walls may severely limit the accuracy and physiological understanding obtained from the simulations.
Therefore, it is essential to utilize CFD solvers with enhanced capabilities, such as moving domain features.

The aim for `OasisMove` was to create a versatile moving domain CFD solver that facilitates open and reproducible
science, with an emphasis on user-friendliness, and geared towards students and researchers. The software extends the
CFD solver `Oasis` [@mortensen2015oasis], which has been thoroughly verified and
validated [@khan2019direct; @bergersen2019fda], and widely used for hemodynamic simulations in the
past [@valen2014mind; @valen2014high; @mancini2019high]. `OasisMove` is implemented in the FEniCS finite element
framework and features an intuitive high-level Python interface with close similarity between the code and the
mathematical equations. Being based on FEniCS, the solver also comes with seamless parallel implementation and
entry-level high-performance computing capabilities.
`OasisMove` fills the aforementioned gaps in moving domain CFD by expressing the Navier-Stokes equations in the
arbitrary Lagrangian-Eulerian formulation, making it suitable for handling moving domains. The solver has undergone
rigorous verification and validation, displaying excellent agreement with theoretical convergence rates, high-resolution
numerical benchmarks, and experimental results [@kjeldsberg2023verified]. As a result, `OasisMove` is a robust and
thoroughly tested solver, which can be a valuable resource for researchers in numerous scientific and engineering
disciplines.

![
In the upper row, the mesh deformation for the 2D vortex problem with moving boundaries is visualized for a half cycle. In the lower row, the instantaneous velocity is shown, where the black arrows indicate the direction of the velocity field, and have been scaled by the velocity magnitude. \label{fig:vortex}
](Figure1.png)

# `OasisMove` in Action

To showcase the capabilities of `OasisMove`, we first examine a Taylor-Green-inspired 2D vortex problem featuring
prescribed oscillating boundaries, which gives rise to a periodic vortex pattern. This problem has an analytical
solution, making it an ideal test for code verification. As time progresses, the boundary undergoes oscillatory motion,
and the inner nodes of the mesh are adjusted accordingly, as depicted in the upper row of \autoref{fig:vortex}. In the
lower row, we present the velocity solution over one half period, with the color indicating the velocity magnitude and
the black arrows indicating the direction of the velocity field. Further analysis and quantitative results from this
test problem can be found in the verification and validation paper [@kjeldsberg2023verified].

![
The z-component of the instantaneous vorticity contours over a full motion cycle for an oscillating cylinder in a free stream. The black axis marks the center of the cylinder at $t = 0$. \label{fig:cylinder}](Figure2.png)

Second, we consider flow around an oscillating cylinder, which is a classic example in CFD used to study the dynamics of
fluid motion around a moving object. In this problem, a cylinder is placed in a free stream, and is prescribed a
sinusoidal vertical motion. The problem is challenging because it involves transient flows, and the oscillation of the
cylinder introduces additional complexity. A robust and accurate moving domain CFD solver must be capable of capturing
the fluid dynamics around the moving object, including the formation and shedding of vortices, as well as the surface
stresses acting on the cylinder. In \autoref{fig:cylinder}, we display snapshots of the flow field over one complete
cycle, illustrating the vortex shedding in the cylinder's wake. Quantitatively, we have previously shown
that `OasisMove` is able to reproduce lift and drag coefficients for this problem with excellent agreement against
reference data [@kjeldsberg2023verified].

![
A visualization of the left atrium velocity field, volume rendered over one cardiac cycle at four time points. We also present the volume curve displayed alongside a black vertical line to indicate the time evolution, shown in the upper right corner of each frame. \label{fig:atrium}
](Figure3.png)

We have previously validated our CFD solver for hemodynamics in a 3D left ventricle, showcasing its accuracy in
simulating intricate flow patterns within the heart [@kjeldsberg2023verified]. To further demonstrate one
of `OasisMove`'s intended applications, we here examine the hemodynamics in a left atrium model, obtained from an
open-source and published dataset [@roneycaroline2022], and pre-processed using VaMPy [@kjeldsberg2023vampy]. CFD
analysis of the left atrium is emerging as a valuable tool for investigating the hemodynamic effects of abnormal cardiac
function, for instance to predict blood clot formation and stroke resulting from atrial fibrillation. \autoref{fig:atrium} portrays the dynamic behavior of the left atrium over a single cardiac cycle at four time points, with the
geometry outline in black, and the volume curve displayed alongside a black vertical line to indicate time progression.
The volume-rendered velocity field demonstrates the inflow through the four pulmonary veins and outflow through the
mitral valve, shedding light on the complex three-dimensional hemodynamics of the left atrium. The figure also
emphasizes the deformation of and flow into the left atrial appendage, indicating that our solver is capable of
generating a plausible flow field under physiological conditions.

# Acknowledgements

This work was supported by the SimCardioTest project (Digital transformation in Health and Care SC1-DTH-06-2020) under
grant agreement No. 101016496 and ERACoSysMed PARIS project under grant agreements No. 643271. The simulations were
performed on the Saga cluster, with resources provided by UNINETT Sigma2 – the National Infrastructure for High
Performance Computing and Data Storage in Norway, grant number nn9249k. We wish to acknowledge the open-source
project [`Oasis`](https://github.com/mikaem/Oasis), and [`FEniCS`](https://fenicsproject.org).

# References
