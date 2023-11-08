# Written by Mikael Mortensen <mikaem@math.uio.no> (2013)
# Edited by Henrik Kjeldsberg <henrik.kjeldsberg@live.no> (2023)

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers


# Override some problem specific parameters
def problem_parameters(NS_parameters, scalar_components, **NS_namespace):
    NS_parameters.update(
        # Mesh parameters
        Nx=50,
        Ny=50,
        # Fluid parameters
        nu=0.001,
        # Simulation parameters
        T=10.0,
        dt=0.005,
        folder="results_driven_cavity",
        # Oasis parameters
        testing=False,
        max_iter=1,
        dynamic_mesh=False,
        save_solution_frequency=5,
        checkpoint=500,
        print_intermediate_info=100,
        velocity_degree=1,
        pressure_degree=1,
        use_krylov_solvers=True)

    scalar_components += ["alfa", "beta"]
    Schmidt["alfa"] = 1.
    Schmidt["beta"] = 10.


# Create a mesh
def mesh(Nx=50, Ny=50, **params):
    m = UnitSquareMesh(Nx, Ny)
    return m


# Specify boundary conditions
def create_bcs(V, **NS_namespace):
    noslip = "std::abs(x[0]*x[1]*(1-x[0]))<1e-8"
    top = "std::abs(x[1]-1) < 1e-8"
    bottom = "std::abs(x[1]) < 1e-8"
    bc0 = DirichletBC(V, 0, noslip)
    bc00 = DirichletBC(V, 1, top)
    bc01 = DirichletBC(V, 0, top)
    bcbeta = DirichletBC(V, 1, bottom)
    return dict(u0=[bc00, bc0],
                u1=[bc01, bc0],
                p=[],
                alfa=[bc00],
                beta=[bcbeta])


def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:
        [bc.apply(x_2[ui]) for bc in bcs[ui]]


def pre_solve_hook(mesh, newfolder, velocity_degree, **NS_namespace):
    # Visualization files
    viz_p, viz_u = get_visualization_writers(newfolder, ['pressure', 'velocity'])

    Vv = VectorFunctionSpace(mesh, 'CG', velocity_degree)
    uv = Function(Vv, name="Velocity")

    return dict(uv=uv, viz_u=viz_u, viz_p=viz_p)


def temporal_hook(viz_u, viz_p, tstep, u_, t, uv, p_, plot_interval, testing, **NS_namespace):
    if tstep % plot_interval == 0 and not testing:
        assign(uv.sub(0), u_[0])
        assign(uv.sub(1), u_[1])

        viz_u.write(uv, t)
        viz_p.write(p_, t)


def theend_hook(u_, uv, mesh, testing, **NS_namespace):
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
