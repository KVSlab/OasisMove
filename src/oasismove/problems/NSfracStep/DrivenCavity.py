# Written by Mikael Mortensen <mikaem@math.uio.no> (2013)
# Edited by Henrik Kjeldsberg <henrik.kjeldsberg@live.no> (2023)

from oasismove.problems.DrivenCavity import *
from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_files


# Override some problem specific parameters
def problem_parameters(NS_parameters, **NS_namespace):
    NS_parameters.update(
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
        save_solution_frequency=1,
        checkpoint=500,
        print_intermediate_info=100,
        velocity_degree=1,
        pressure_degree=1,
        use_krylov_solvers=True,
        max_error=1e-8)


# Specify boundary conditions
def create_bcs(V, **NS_namespace):
    bc0 = DirichletBC(V, 0, noslip)
    bc00 = DirichletBC(V, 1, top)
    bc01 = DirichletBC(V, 0, top)
    return dict(u0=[bc00, bc0],
                u1=[bc01, bc0],
                p=[])


def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:
        [bc.apply(x_2[ui]) for bc in bcs[ui]]


def pre_solve_hook(mesh, newfolder, velocity_degree, **NS_namespace):
    # Visualization files
    viz_p, viz_u = get_visualization_files(newfolder)

    Vv = VectorFunctionSpace(mesh, 'CG', velocity_degree)
    return dict(uv=Function(Vv), viz_u=viz_u, viz_p=viz_p)


def temporal_hook(viz_u, viz_p, tstep, u_, t, uv, p_, plot_interval, testing, **NS_namespace):
    if tstep % plot_interval == 0 and not testing:
        assign(uv.sub(0), u_[0])
        assign(uv.sub(1), u_[1])

        viz_u.write(uv, t)
        viz_p.write(p_, t)


def theend_hook(u_, p_, tstep, save_solution_frequency, uv, mesh, testing, **NS_namespace):
    u_norm = norm(u_[0].vector())
    if MPI.rank(MPI.comm_world) == 0 and testing:
        print("Velocity norm = {0:2.6e}".format(u_norm))

    if not testing:
        try:
            from fenicstools import StreamFunction
            psi = StreamFunction(uv, [], mesh, use_strong_bc=True)
            plot(psi, title='Streamfunction')
            import matplotlib.pyplot as plt
            plt.show()
        except ImportError:
            pass
