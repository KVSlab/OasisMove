# Parameters

Every problem will contain parameters that are both specific to the problem, but also general to Oasis. All parameters
will be collected in the `NS_parameters` dictionary, and accessible in every OasisMove [hook](gs:hooks). You can also
adjust the parameters directly though the command line, by passing them as command line arguments.

## Overview of parameters

To get an impression of what parameters you might come accross in the problem files included in OasisMove, we will
consider the simplest rigid flow problem as an example: the
default [DrivenCavity](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/DrivenCavity.py)
problem. By inspecting the problem file, the parameters are defined in the start of the script:

```python
NS_parameters.update(
    # Mesh parameters
    Nx=50,
    Ny=50,

    # Fluid parameters
    nu=0.001,

    # Simulation parameters
    T=1.0,
    dt=0.005,
    folder="results_driven_cavity",

    # Oasis parameters
    testing=False,
    max_iter=2,
    dynamic_mesh=False,
    save_solution_frequency=5,
    checkpoint=500,
    print_intermediate_info=100,
    velocity_degree=1,
    pressure_degree=1,
    use_krylov_solvers=True
)
```

In the OasisMove problem files, we distinguish between `Mesh`, `Fluid`, `Simulation`, and `Oasis` parameters. In
addition, some problems have parameters that are specific for the particular problems. The full dictionary of default
parameters can be found in the top of the
in [problems/NSfracStep/\_\_init\_\_.py](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/__init__.py)
script.

### `Mesh` parameters

Mesh parameters are related to the computational mesh, and typically the resolution. In regard to the `DrivenCavity`
problem, the mesh parameters are:

```python
# Mesh parameters
Nx = 50
Ny = 50
```

By adjusting `Nx` and `Ny` we adjust the resolution of the computational domain, which for the `DrivenCavity` problem is
a 2D unit square mesh consisting of triangular elements. The parameters `Nx` and `Ny` determine the number of uniform
cell intervals in the $x$ and $y$ direction, respectively. These parameters can also be used to compute the totan number
of triangular cells in 2D: or tetrehedral cells (in 3D) the computational domain consists of:

```{math}
    \text{Number of cells (2D) } &= 2 \times N_x \times N_y             \\
    \text{Number of cells (3D) } &= 6 \times N_x \times N_y  \times N_z 
```
### `Fluid` parameters

Fluid parameters are related to the fluid properties. In regard to the `DrivenCavity`
problem, the only fluid parameter is the kinematic viscosity ($\nu$):

```python
# Fluid parameters
nu = 0.001
```

Changing the kinematic viscosity will directly influence the Reynolds number, the ratio between inertial and viscous
forces defined as:

```{math}
:label: eq-reynolds
\begin{align}
    \text{Re } = \frac{u L}{\nu}, 
\end{align}
```

where $u$ is the flow speed and $L$ is a characteristic linear dimension.

### `Simulation` parameters

Simulation parameters are related to the computational fluid dynamics simulation. In regard to the `DrivenCavity`
problem, the simulation parameters are:

```python
# Simulation parameters
T = 1.0
dt = 0.005
folder = "results_driven_cavity"
```

Here, `T` is the total simulation time, usually in seconds or milliseconds depending on the application and scale of the
problem. The number of time steps is controlled by the value of `T`, and the time step size `dt`, which are related as
follows:

```{math}
:label: eq-dt
\begin{align}
\text{Number of time steps } = \frac{T}{\Delta t}.
\end{align}
```

Hence, for the default `DrivenCavity` problem, the simulation is run for `T` / `dt` = 200 time steps. The `folder`
parameter determines where the results are saved for the particular problem, and is intentionally different for every
problem to easier distinguish the results. Also, the `**NS_namespace` parameter `newfolder` is a derivative of
the `folder` location:

```console
newfolder = folder/data/[RUN NUMBER]/ 
```

and is used to determine where the solution files, which are stored to:

```console
newfolder/Solutions 
```

### `Oasis` parameters

Oasis parameters are miscellaneous parameters related to solver performance, printing simulation info, or parameters
inherited from classical Oasis. In regard to the `DrivenCavity` problem, the Oasis parameters are:

```python
# Oasis parameters
testing = False
max_iter = 2
dynamic_mesh = False
save_solution_frequency = 5
checkpoint = 500
print_intermediate_info = 100
velocity_degree = 1
pressure_degree = 1
use_krylov_solvers = True
```

The `testing` flag is there to notify the program where or not it is run as a test. The `max_iter` parameter is set to 2
for all moving domain simulations in order to ensure convergence of velocity and pressure, based on empirical evidence.
By default `dynamic_mesh` is `False`, and controlls when the mesh equations are solved followed by mesh coordinate
deformation. The `save_solution_frequency` parameter determines how often the velocity and pressure solution are stored.
With the proposed value of `5`, the simulation will store the solution every 5th time step. `print_intermediate_info`
determines the frequency of printing simulation info, here set to every 100 time step. Furthermore,
the `velocity_degree` and `pressure_degree` parameters determine the finite element polynomial approximation order.
If `use_krylov_solvers` is `True` the solver will use the PETSc Krylov solver, instead of a standard LU solver.

## Changing parameters

Parameters can either be changed directly in the Python problem file, or directly though the command line, by passing
them as command line arguments. For instance, if we want to increase the mesh resolution and extend the simulation end
time to $T=20$ for the `DrivenCavity` problem, you can execute the following command:

``` console
$ oasism NSfracStepMove problem=DrivenCavity Nx=100 Ny=100 T=20  
```

Similarly, you can change any of parameters mention here, the ones listed in
the [problems/NSfracStep/\_\_init\_\_.py](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/__init__.py)
script, or that are specified in the Python problem script.
