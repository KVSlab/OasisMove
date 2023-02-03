# Hooks
This is an updated version of the Oasis wiki's [Hooks](https://github.com/mikaem/Oasis/wiki/Hooks), adjusted with the new functionality in OasisMove.
The Navier Stokes solvers `NSfracStep` or `NSCoupled` both need to import a mesh, some parameters and a number of solver specific or user defined hooks from the solvers and problems submodules respectively. The problem hooks are called at certain strategic locations in the solvers. All problem hooks have default versions in [problems/NSfracStep/\_\_init\_\_.py](https://github.com/mikaem/Oasis/blob/master/problems/NSfracStep/__init__.py) and [problems/NSCoupled/\_\_init\_\_.py](https://github.com/mikaem/Oasis/blob/master/problems/NSCoupled/__init__.py). The default hooks for `NSfracStep` with brief descriptions are basically:

```python
def body_force(mesh, **NS_namespace):
    """Specify body force"""
    return Constant((0,)*mesh.geometry().dim())

def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    return dict((ci, Constant(0)) for ci in scalar_components)
    
def initialize(**NS_namespace):
    """Initialize solution."""
    pass

def create_bcs(sys_comp, **NS_namespace):
    """Return dictionary of Dirichlet boundary conditions."""
    return dict((ui, []) for ui in sys_comp)

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
All hooks may be overloaded in the problem module because the problem imports from the `__init__.py` file initially.

The solver also have a list of functions that may be overloaded. The default functions required by the `NSfracStep` solver is

```python
__all__ = ["assemble_first_inner_iter", "velocity_tentative_assemble",
           "velocity_tentative_solve", "pressure_assemble", 
           "pressure_solve", "velocity_update", "scalar_assemble", 
           "scalar_solve", "get_solvers", "setup", 
           "print_velocity_pressure_info"]

def setup(**NS_namespace):
    """Set up all equations to be solved."""
    return {}

def get_solvers(**NS_namespace):
    """Return 4 linear solvers. 
    
    We are solving for
       - tentative velocity
       - pressure correction
       - velocity update (unless lumping is switched on)
       
       and possibly:       
       - scalars
            
    """        
    return (None, )*4

def assemble_first_inner_iter(**NS_namespace):
    """Called first thing on a new velocity/pressure iteration."""
    pass

def velocity_tentative_solve(**NS_namespace):
    """Linear algebra solve of tentative velocity component."""
    pass

def velocity_tentative_assemble(**NS_namespace):
    """Assemble remaining system for tentative velocity component."""
    pass
    
def pressure_assemble(**NS_namespace):
    """Assemble rhs of pressure equation."""
    pass

def pressure_solve(**NS_namespace):
    """Solve pressure equation."""    
    pass

def velocity_update(**NS_namespace):
    """Update the velocity after finishing pressure velocity iterations."""
    pass

def scalar_assemble(**NS_namespace):
    """Assemble scalar equation."""
    pass

def scalar_solve(**NS_namespace):
    """Solve scalar equation."""
    pass

```

Note the special keyword argument `**NS_namespace`. This is actually the entire namespace of the `NSfracStep` solver. The hooks are called using the following notation:
```python
initialize(**vars())

temporal_hook(**vars())

```
That is, the solver sends its entire namespace, vars(), to the problem as keyword arguments. In the problem module you may now access anything from the Navier Stokes solver when initializing. The hook can be used to set the velocity and pressure fields through for example the following code in the problem module:
```python
init_fields = {"u0": 0, "u1": 1, "u2": 2, "p": 0}
def initialize(q_, q_1, q_2, **NS_namespace):
    for ui in q_:
        q_ [ui].vector()[:] = init_fields[ui]
    for ui in q_1:
        q_1[ui].vector()[:] = init_fields[ui]
        q_2[ui].vector()[:] = init_fields[ui]
```
Note that `q_` and `q_1` are dictionaries declared in the solver that hold the solution at current and previous timesteps. For example, for 2D problems `q_` is declared as
```python
V = FunctionSpace(mesh, "CG", velocity_degree)
Q = FunctionSpace(mesh, "CG", pressure_degree)
q_ = {"u0": Function(V), "u1": Function(V), "p": Function(Q)}
```
and similarily for `q_1` and `q_2`, only that these do not contain any pressures, since the pressure is only stored for one timestep. The dictionaries are accessed through arguments of the `initialize` hook. Any variable that you want to use is simply written out, the rest is soaked up by `**NS_namespace`.

