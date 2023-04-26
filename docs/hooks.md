(gs:hooks)=

# Hooks

This is an updated version of the Oasis wiki's [Hooks](https://github.com/mikaem/Oasis/wiki/Hooks), adjusted with the
new functionality of OasisMove. The Navier-Stokes solvers `NSfracStep`, `NSCoupled`, and `NSfracStepMove` need to import
a mesh, some parameters and a number of solver specific or user defined hooks from the solvers and problems submodules
respectively. The problem hooks are called at certain strategic locations in the solvers. All problem hooks have default
versions
in [problems/NSfracStep/\_\_init\_\_.py](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/__init__.py)
and [problems/NSCoupled/\_\_init\_\_.py](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSCoupled/__init__.py)
.

The default hooks for `NSfracStepMove` with brief descriptions are:

```python
def body_force(mesh, **NS_namespace):
    """Specify body force"""
    return Constant((0,) * mesh.geometry().dim())


def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    return dict((ci, Constant(0)) for ci in scalar_components)


def initialize(**NS_namespace):
    """Initialize solution."""
    pass


def pre_boundary_condition(**NS_namespace):
    """Called after defining the function spaces and functions,
       but before the boundary conditions are created.
    """
    return dict(boundary=None)


def create_bcs(sys_comp, **NS_namespace):
    """Return dictionary of Dirichlet boundary conditions."""
    return dict((ui, []) for ui in sys_comp)


def update_boundary_conditions(**NS_namespace):
    """Called just prior of assembling mesh and/or fluid equations.
    """


def mesh_velocity_hook(**NS_namespace):
    """Called just prior to solving for mesh velocity."""
    pass


def velocity_tentative_hook(**NS_namespace):
    """Called just prior to solving for tentative velocity."""
    pass


def pressure_hook(**NS_namespace):
    """Called prior to pressure solve."""
    pass


def scalar_hook(**NS_namespace):
    """Called prior to scalar solve."""
    pass


def start_timestep_hook(**NS_parameters):
    """Called at start of new timestep"""
    pass


def temporal_hook(**NS_namespace):
    """Called at end of a timestep."""
    pass


def pre_solve_hook(**NS_namespace):
    """Called just prior to entering time-loop. Must return a dictionary."""
    return {}


def theend_hook(**NS_namespace):
    """Called at the very end."""
    pass

```

All hooks may be overloaded in the problem module because the problem imports from the `__init__.py` file initially. For
information of the classical Oasis hooks, and list of solver hooks we refer to the
Oasis [wiki](https://github.com/mikaem/Oasis/wiki/Hooks).

## Understanding `NS_namespace`

Throughout the Python scripts you will notice the special keyword argument `**NS_namespace`. This is actually the entire
namespace of the `NSfracStepMove` solver, where the hooks are called using the following notation:

```python
initialize(**vars())

temporal_hook(**vars())
```

That is, the solver sends its entire namespace, vars(), to the problem as keyword arguments. In the problem module you
may now access anything from the Navier-Stokes solver when initializing. The hook can be used to set the velocity and
pressure fields through for example the following code in the problem module, which has been done in
the [MovingVortex](https://github.com/KVSlab/OasisMove/blob/main/src/oasismove/problems/NSfracStep/MovingVortex.py)
problem:

```python
init_fields = {"u0": 0, "u1": 1, "u2": 2, "p": 0}


def initialize(q_, q_1, q_2, **NS_namespace):
    for ui in q_:
        q_[ui].vector()[:] = init_fields[ui]
    for ui in q_1:
        q_1[ui].vector()[:] = init_fields[ui]
        q_2[ui].vector()[:] = init_fields[ui]
```

Note that `q_` and `q_1` are dictionaries declared in the solver that hold the solution at current and previous
timesteps. For example, for 2D problems `q_` is declared as

```python
V = FunctionSpace(mesh, "CG", velocity_degree)
Q = FunctionSpace(mesh, "CG", pressure_degree)
q_ = {"u0": Function(V), "u1": Function(V), "p": Function(Q)}
```

and similarily for `q_1` and `q_2`, only that these do not contain any pressures, since the pressure is only stored for
one timestep. The dictionaries are accessed through arguments of the `initialize` hook. Any variable that you want to
use is simply written out, the rest is soaked up by `**NS_namespace`.

