# OasisMove - A verified and validated Python/FEniCS-based CFD solver for moving domains 

[OasisMove](https://github.com/KVSlab/OasisMove) is a high-level/high-performance open-source Navier-Stokes solver written in Python/[FEniCS](https://fenicsproject.org/), and is an extension of the computational fluid dynamics (CFD) solver [Oasis](https://github.com/mikaem/Oasis).
In OasisMove the Navier-Stokes equations are expressed in the arbitrary Lagrangian-Eulerian formulation, which is suitable for handling moving domains.
Through verification the solver has shown to follow theoretical convergence rates, begin second order accurate in time, and second and third order accurate in space with P1/P1 and P2/P1 finite elements.
Through validation the solver has showed good agreement with existing benchmark results, and demonstrates the ability of the solver to capture vortex patterns in transitional and turbulent-like flow regimes.

## Contents
```{tableofcontents}
```
