import os
import pickle

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers


# Override some problem specific parameters
def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    """
    Problem file for running CFD simulation for a stenosis model with a slight eccentricity in a 2D domain.
    The stenotic and eccentricity model is inspired by the 3D stenotic flow research by  Varghese et al.[1], where we
    have reduced the length of the stenosis considerably on the right-hand side. The stenosis model is prescribed a
    parabolic inlet profile with maximum velocity U0, and an open outlet with a zero pressure condition.
    The problem is intended to demonstrate backflow stabilization proposed by Moghadam et al. [2], which is required
    with the default paramaters, run at Re=2540 and controlled by the parameters' backflow_facet and backflow_beta,
    representing a list of boundaries to apply the stabilization to, and the stabilization strength, respectively.

    [1] Varghese, S. S., Frankel, S. H., & Fischer, P. F. (2007). Direct numerical simulation of stenotic flows. Part 1.
    Steady flow. Journal of Fluid Mechanics, 582, 253-280.
    [2] Esmaily Moghadam, M., Bazilevs, Y., Hsia, T. Y., Vignon-Clementel, I. E., Marsden, A. L., & Modeling of
    Congenital Hearts Alliance (MOCHA). (2011). A comparison of outlet boundary treatments for prevention of backflow
    divergence with relevance to blood flow simulations. Computational Mechanics, 48, 277-291.
    """
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(os.getcwd(), restart_folder)
        f = open(path.join(path.dirname(path.abspath(__file__)), restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)
    else:
        NS_parameters.update(
            #  Backflow parameters
            backflow_facets=[3],  # Outlet with id=3 needs backflow stabilization.
            backflow_beta=0.2,  # Strength of backflow stabilization
            # Problem specific parameters
            x0=-3,  # Min x-coodrinate of the domain
            x1=3,  # Max x-coodrinate of the domain
            N=50,  # Mesh resolution
            D=6.35,  # Diameter
            U0=4,  # Maximum inlet velocity
            nu=0.01,  # Kinetmatic viscosity
            # Simulation parameters
            T=15,  # Time
            dt=1e-2,  # Time step size
            folder="results_stenosis",
            # Oasis paramters
            max_iter=1,
            dynamic_mesh=False,
            save_solution_frequency=5,
            checkpoint=500,
            print_intermediate_info=100,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True
        )


def mesh(N, D, x0, x1, dt, U0, **params):
    # Define horizontal and verticl resolution
    Ny = N
    Nx = (x1 - x0) * Ny

    # Stenosis function
    def S(x, s0=0.25):
        L = 2 * D
        return 1 / 2 * D * (1 - s0 * (1 + np.cos(2 * np.pi * x / L)))

    # Eccentricity function
    def E(x, s0=0.25):
        L = 2 * D
        return 0.1 * D * s0 * (1 + np.cos(2 * np.pi * x / L))

    # Create stenosis mesh
    x0 = x0 * D
    y0 = -D / 2

    x1 = x1 * D
    y1 = D / 2
    mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), Nx, Ny)

    # Squeeze towards walls
    c = mesh.coordinates()

    x = c[:, 0]
    y = c[:, 1]

    for i in np.arange(c.shape[0]):
        if -D <= x[i] <= D:
            y[i] = S(x[i]) * y[i] / (D / 2) + E(x[i])

    print_mesh_information(mesh, dt, U0, dim=2)
    return mesh


def pre_boundary_condition(mesh, D, x0, x1, **NS_namespace):
    # Define boundaries
    x0 = x0 * D
    x1 = x1 * D
    wall = AutoSubDomain(
        lambda x, b: (b and not near(x[0], x0 - 1000 * DOLFIN_EPS) and not near(x[0], x1 + 1000 * DOLFIN_EPS)))
    inlet = AutoSubDomain(lambda x, b: (b and (near(x[0], x0))))
    outlet = AutoSubDomain(lambda x, b: (b and (near(x[0], x1))))

    # Mark geometry boundaries
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    wall.mark(boundary, 1)
    inlet.mark(boundary, 2)
    outlet.mark(boundary, 3)

    return dict(boundary=boundary)


def create_bcs(V, Q, D, U0, boundary, sys_comp, **NS_namespace):
    info_red("Creating boundary conditions")

    # Walls
    bc_wall = DirichletBC(V, Constant(0), boundary, 1)

    # Inlet
    parabolic = Expression("U0*(1 - 4 * x[1]*x[1] / (D*D) )", U0=U0, D=D, degree=2)
    bc_in_x = DirichletBC(V, parabolic, boundary, 2)
    bc_in_y = DirichletBC(V, Constant(0), boundary, 2)

    # Outlet
    bc_out = DirichletBC(Q, Constant(0), boundary, 3)

    bcs = dict((ui, []) for ui in sys_comp)
    bcs['u0'] = [bc_wall, bc_in_x]
    bcs['u1'] = [bc_wall, bc_in_y]
    bcs['p'] = [bc_out]

    return bcs


def pre_solve_hook(V, mesh, newfolder, velocity_degree, **NS_namespace):
    # Get visualization files
    viz_p, viz_u = get_visualization_writers(newfolder, ['pressure', 'velocity'])
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(Vv, name="u")

    return dict(viz_u=viz_u, viz_p=viz_p, u_vec=u_vec)


def temporal_hook(tstep, u_vec, u_, save_solution_frequency, viz_p, viz_u, p_, t, T, **NS_namespace):
    if tstep % save_solution_frequency == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])

        viz_u.write(u_vec, t)
        viz_p.write(p_, t)

    if tstep % 10 == 0:
        u_x = u_[0].vector().get_local()
        max_vel = max(u_x)
        mean_vel = np.mean(u_x)
        info_green(
            'Time = {0:2.4e}, timestep = {1:6d}, max velocity={2:2.2f} mean velocity={3:2.3f} End time = {4:2.4e}'
            .format(t, tstep, max_vel, mean_vel, T)
        )
