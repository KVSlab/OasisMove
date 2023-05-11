#!/usr/bin/env python

__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2013-11-06"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

"""
This module implements a generic form of the fractional step method for
solving the incompressible Navier-Stokes equations. There are several
possible implementations of the pressure correction and the more low-level
details are chosen at run-time and imported from any one of:

  solvers/NSfracStep/IPCS_ABCN.py    # Implicit convection
  solvers/NSfracStep/IPCS_ABE.py     # Explicit convection
  solvers/NSfracStep/IPCS.py         # Naive implict convection
  solvers/NSfracStep/BDFPC.py        # Naive Backwards Differencing IPCS in rotational form
  solvers/NSfracStep/BDFPC_Fast.py   # Fast Backwards Differencing IPCS in rotational form
  solvers/NSfracStep/Chorin.py       # Naive

The naive solvers are very simple and not optimized. They are intended
for validation of the other optimized versions. The fractional step method
can be used both non-iteratively or with iterations over the pressure-
velocity system.

The velocity vector is segregated, and we use three (in 3D) scalar
velocity components.

Each new problem needs to implement a new problem module to be placed in
the problems/NSfracStep folder. From the problems module one needs to import
a mesh and a control dictionary called NS_parameters. See
problems/NSfracStep/__init__.py for all possible parameters.

"""
import importlib
import pickle
from pprint import pprint

from oasismove.common import *

commandline_kwargs = parse_command_line()

# Find the problem module
default_problem = 'DrivenCavity'
problemname = commandline_kwargs.get('problem', default_problem)
problemspec = importlib.util.find_spec('.'.join(('oasismove.problems.NSfracStep', problemname)))
if problemspec is None:
    problemspec = importlib.util.find_spec(problemname)
if problemspec is None:
    raise RuntimeError(problemname + ' not found')

# Import the problem module
print('Importing problem module ' + problemname + ':\n' + problemspec.origin)
problemmod = importlib.util.module_from_spec(problemspec)
problemspec.loader.exec_module(problemmod)

vars().update(**vars(problemmod))

# Update problem spesific parameters
problem_parameters(**vars())

if MPI.rank(MPI.comm_world) == 0:
    print("=== Starting simulation for {} problem ===".format(problemname))
    print("Running with the following parameters:")
    pprint(NS_parameters)

# Update current namespace with NS_parameters and commandline_kwargs ++
vars().update(post_import_problem(**vars()))

# Use t and tstep from stored paramteres if restarting
previous_velocity_degree = velocity_degree
if restart_folder is not None:
    f = open(path.join(path.abspath(restart_folder), 'params.dat'), 'rb')
    params = pickle.load(f)
    f.close()
    t = params["t"]
    tstep = params["tstep"]
    previous_velocity_degree = params["velocity_degree"]

# Import chosen functionality from solvers
solver = importlib.import_module('.'.join(('oasismove.solvers.NSfracStep', solver)))
vars().update({name: solver.__dict__[name] for name in solver.__all__})

# Create lists of components solved for
dim = mesh.geometry().dim()
u_components = ['u' + str(x) for x in range(dim)]
sys_comp = u_components + ['p'] + scalar_components
uc_comp = u_components + scalar_components

# Set up initial folders for storing results
newfolder, tstepfiles = create_initial_folders(**vars())

# Declare FunctionSpaces and arguments
V = Q = FunctionSpace(mesh, 'CG', velocity_degree,
                      constrained_domain=constrained_domain)
if velocity_degree != pressure_degree:
    Q = FunctionSpace(mesh, 'CG', pressure_degree,
                      constrained_domain=constrained_domain)

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)
u_w = TrialFunction(V)
v_w = TestFunction(V)

# Use dictionary to hold all FunctionSpaces
VV = dict((ui, V) for ui in uc_comp)
VV['p'] = Q

# Create dictionaries for the solutions at three timesteps
q_ = dict((ui, Function(VV[ui], name=ui)) for ui in sys_comp)
q_1 = dict((ui, Function(VV[ui], name=ui + "_1")) for ui in sys_comp)
q_2 = dict((ui, Function(V, name=ui + "_2")) for ui in u_components)

# Read in previous solution if restarting
init_from_restart(**vars())

# Create dictionary for the wall solution at two timestep
w_ = dict((ui, Function(V, name=ui)) for ui in u_components)

# Create vectors of the segregated velocity components
u_ = as_vector([q_[ui] for ui in u_components])  # Velocity vector at t
u_1 = as_vector([q_1[ui] for ui in u_components])  # Velocity vector at t - dt
u_2 = as_vector([q_2[ui] for ui in u_components])  # Velocity vector at t - 2*dt

# Create vectors of the segregated mesh velocity compinents
wu_ = as_vector([w_[ui] for ui in u_components])  # Mesh velocity at t

# Adams Bashforth projection of velocity at t - dt/2
U_AB = 1.5 * u_1 - 0.5 * u_2

# Create short forms for accessing the solution vectors
x_ = dict((ui, q_[ui].vector()) for ui in sys_comp)  # Solution vectors t
x_1 = dict((ui, q_1[ui].vector()) for ui in sys_comp)  # Solution vectors t - dt
x_2 = dict((ui, q_2[ui].vector()) for ui in u_components)  # Solution vectors t - 2*dt
wx_ = dict((ui, w_[ui].vector()) for ui in u_components)

# Create vectors to hold rhs of equations
b = dict((ui, Vector(x_[ui])) for ui in sys_comp)  # rhs vectors (final)
b_tmp = dict((ui, Vector(x_[ui])) for ui in sys_comp)  # rhs temp storage vectors
bw = dict((ui, Vector(wx_[ui])) for ui in u_components)  # rhs vectors (final)
bw_tmp = dict((ui, Vector(wx_[ui])) for ui in u_components)  # rhs temp storage vectors

# Short forms pressure and scalars
p_ = q_['p']  # pressure at t
p_1 = q_1['p']  # pressure at t - dt
dp_ = Function(Q)  # pressure correction
for ci in scalar_components:
    exec("{}_   = q_ ['{}']".format(ci, ci))
    exec("{}_1  = q_1['{}']".format(ci, ci))

print_solve_info = use_krylov_solvers and krylov_solvers['monitor_convergence']

# Anything problem specific
vars().update(pre_boundary_condition(**vars()))

# Boundary conditions
bcs = create_bcs(**vars())

# LES setup
# exec("from oasis.solvers.NSfracStep.LES.{} import *".format(les_model))
lesmodel = importlib.import_module('.'.join(('oasismove.solvers.NSfracStep.LES', les_model)))
vars().update({name: lesmodel.__dict__[name] for name in lesmodel.__all__})

vars().update(les_setup(**vars()))

# Non-Newtonian setup
# exec("from oasis.solvers.NSfracStep.NNModel.{} import *".format(nn_model))
nnmodel = importlib.import_module('.'.join(('oasismove.solvers.NSfracStep.NNModel', nn_model)))
vars().update({name: nnmodel.__dict__[name] for name in nnmodel.__all__})

vars().update(nn_setup(**vars()))

# Initialize solution
initialize(**vars())

#  Fetch linear algebra solvers
u_sol, p_sol, w_sol, c_sol = get_solvers(**vars())

# Get constant body forces
f = body_force(**vars())
assert (isinstance(f, Coefficient))
b0 = dict((ui, assemble(v * f[i] * dx)) for i, ui in enumerate(u_components))
bw0 = dict((ui, assemble(v * Constant(0) * dx)) for i, ui in enumerate(u_components))

# Get scalar sources
fs = scalar_source(**vars())
for ci in scalar_components:
    assert (isinstance(fs[ci], Coefficient))
    b0[ci] = assemble(v * fs[ci] * dx)

# Preassemble and allocate
vars().update(setup(**vars()))

# Anything problem specific
vars().update(pre_solve_hook(**vars()))

tx = OasisTimer('Timestep timer')
tx.start()
stop = False
total_timer = OasisTimer("Start simulations", print_solve_info)

max_tstep = 10 #if restart_folder is None else tstep + 10

print("Saving results to: {}".format(newfolder))
while t < (T - tstep * DOLFIN_EPS) and not stop:
    t_velocity = 0
    t += dt
    tstep += 1
    inner_iter = 0
    udiff = array([1e8])  # Norm of velocity change over last inner iter
    num_iter = max(iters_on_first_timestep, max_iter) if tstep <= max_tstep else max_iter

    start_timestep_hook(**vars())

    update_boundary_conditions(**vars())

    # Solve for mesh velocity and move mesh with prescribed motion
    if dynamic_mesh:
        t0 = OasisTimer("Mesh equations")
        for i, ui in enumerate(u_components):
            mesh_velocity_assemble(**vars())
            mesh_velocity_hook(**vars())
            mesh_velocity_solve(**vars())
        t0.stop()

        b0 = dict((ui, assemble(v * f[i] * dx)) for i, ui in enumerate(u_components))
        A_cache.update_t(t)

    while udiff[0] > max_error and inner_iter < num_iter and compute_velocity_and_pressure:
        inner_iter += 1

        if inner_iter == 1:
            les_update(**vars())
            nn_update(**vars())
            assemble_first_inner_iter(**vars())

        udiff[0] = 0.0
        t0 = OasisTimer("Tentative velocity")
        for i, ui in enumerate(u_components):
            t1 = OasisTimer('Solving tentative velocity ' + ui, print_solve_info)
            velocity_tentative_assemble(**vars())
            velocity_tentative_hook(**vars())
            velocity_tentative_solve(**vars())
            t1.stop()
        t0.stop()
        if print_solve_info:
            info_red("Tentative velocity (Total) Elapsed time: {:.5f}".format(t0.elapsed()[0]))

        t0 = OasisTimer("Pressure solve", print_solve_info)
        pressure_assemble(**vars())
        pressure_hook(**vars())
        pressure_solve(**vars())
        t0.stop()
        if print_solve_info:
            info_red("Pressure (Total) Elapsed time: {:.5f}".format(t0.elapsed()[0]))

        print_velocity_pressure_info(**vars())
    # Update velocity
    if compute_velocity_and_pressure:
        t0 = OasisTimer("Velocity update", print_solve_info)
        velocity_update(**vars())
        t0.stop()
        if print_solve_info:
            info_red("Elapsed time: {:.5f}".format(t0.elapsed()[0]))

        # Solve for scalars
        if len(scalar_components) > 0:
            scalar_assemble(**vars())
            for ci in scalar_components:
                t1 = OasisTimer('Solving scalar {}'.format(ci), print_solve_info)
                scalar_hook(**vars())
                scalar_solve(**vars())
                t1.stop()

    temporal_hook(**vars())

    # Save solution if required and check for killoasis file
    stop = save_solution(**vars())

    if compute_velocity_and_pressure:
        # Update to a new timestep
        for ui in u_components:
            x_2[ui].zero()
            x_2[ui].axpy(1.0, x_1[ui])
            x_1[ui].zero()
            x_1[ui].axpy(1.0, x_[ui])

        for ci in scalar_components:
            x_1[ci].zero()
            x_1[ci].axpy(1., x_[ci])

    # Print some information
    if tstep % print_intermediate_info == 0:
        toc = tx.stop()
        info_green('Time = {0:2.4e}, timestep = {1:6d}, End time = {2:2.4e}'.format(t, tstep, T))
        info_red('Total computing time on previous {0:d} timesteps = {1:f}'.format(
            print_intermediate_info, toc))
        list_timings(TimingClear.clear, [TimingType.wall])
        tx.start()

    # AB projection for pressure on next timestep
    if AB_projection_pressure and t < (T - tstep * DOLFIN_EPS) and not stop:
        x_['p'].axpy(0.5, dp_.vector())

total_timer.stop()
list_timings(TimingClear.keep, [TimingType.wall])
info_red('Total computing time = {0:f}'.format(total_timer.elapsed()[0]))
oasis_memory('Final memory use ')
total_initial_dolfin_memory = MPI.sum(MPI.comm_world, initial_memory_use)
info_red('Memory use for importing dolfin = {} MB (RSS)'.format(
    total_initial_dolfin_memory))
info_red('Total memory use of solver = ' +
         str(oasis_memory.memory - total_initial_dolfin_memory) + " MB (RSS)")

if restart_folder is not None:
    merge_visualization_files(**vars())

# Final hook
theend_hook(**vars())
