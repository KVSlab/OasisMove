(tutorial:interpolate)=

# Interpolating to a higher function space

OasisMove provides the ability to interpolate a solution to a higher function space. For instance, interpolating a
$\mathbb{P_1}$ solution to $\mathbb{P_2}$ function space. This is particularly useful when you want to increase the
resolution of your solution, or when you need a higher degree of accuracy. Note that this feature is currently
implemented for the velocity field only.

## Interpolating the velocity to a higher function space

The interpolation technique in OasisMove is closely related to [restarting a simulation](tutorial:restart). To
demonstrate interpolation of the velocity from $\mathbb{P_1}$ to $\mathbb{P_2}$ function space, we will consider
the [`DrivenCavity`](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/DrivenCavity.py)
problem. By default, the problem runs for $T=10$ using $\mathbb{P_1}$ elements, determined by the
argument `velocity_degree=1`. Assume we have run the problem using the default parameters:

``` console
$ oasismove NSfracStep solver=IPCS_ABCN problem=DrivenCavity 
```

To restart and continue the simulation until $T=20$ using $\mathbb{P_2}$ elements, we may simply pass
the `velocity_degree=2` argument, in addition to the required parameters for restarting a simulation:

``` console
$ oasismove NSfracStep solver=IPCS_ABCN problem=DrivenCavity velocity_degree=2 restart_folder=results_driven_cavity/data/1/Checkpoint T=20 
```

Note that the ability to interpolate to a higher function space can be computation-intensive, especially for large
simulations and higher function spaces.
