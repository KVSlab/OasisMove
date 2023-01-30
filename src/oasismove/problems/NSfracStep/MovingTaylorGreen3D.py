from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers


def problem_parameters(NS_parameters, NS_expressions, **NS_namespace):
    """
    Problem file for running CFD simulation for the moving Taylor-Green problem in 3D as described by Fehn et al.[1].
    The problem solves the N-S equations in the absence of body forces, and is commonly used to study transitional and
    turbulent flows. The problem initializes the solution at the two previous time steps, and applies periodic boundary
    condition on the domain walls in all coordinate directions. The domain boundaries are moving for the given mesh
    motion, but the mesh deformation is defined periodically in order to ensure consistency with the periodic boundary
    conditions. The main simulation parameters are the ampliture A, and period time duration T_G.

    [1] Fehn, N., Heinz, J., Wall, W. A., & Kronbichler, M. (2021). High-order arbitrary Lagrangian–Eulerian
    discontinuousGalerkin methods for the incompressible Navier–Stokes equations.
    Journal of Computational Physics, 430, 110040.
    """

    # Override some problem specific parameters
    T_G = 20
    L = 2 * np.pi
    A0 = np.pi / 6
    Re = 1600
    nu = 1 / Re
    dt = 0.005
    T = 5
    recursive_update(NS_parameters, dict(
        # Mesh params
        L=L,
        A0=A0,
        T_G=T_G,
        # Problem params
        nu=nu,
        T=T,
        dt=dt,
        Nx=32,
        Ny=32,
        Nz=32,
        folder="results_moving_taylor_green_3d",
        max_iter=2,
        velocity_degree=1,
        save_step=20E10,
        save_solution_frequency=5,
        dynamic_mesh=True,
        checkpoint=10000,
        plot_interval=100000,
        use_krylov_solvers=True))

    NS_expressions.update(dict(
        constrained_domain=PeriodicDomain(),
        kin=zeros(1),
        initial_fields_w=dict(
            w0="2 * pi * A / T_G * cos(2 * pi * t / T_G) * sin(2 * pi * (x[1] + L/2)/L) * sin(2*pi* (x[2] + L/2)/L)",
            w1="2*pi*A/T_G*cos(2*pi*t/T_G) * sin(2*pi * (x[0] + L/2)/L) * sin(2*pi* (x[2] + L/2)/L)",
            w2="2*pi*A/T_G*cos(2*pi*t/T_G) * sin(2*pi * (x[0] + L/2)/L) * sin(2*pi* (x[1] + L/2)/L)"
        ),
        initial_fields=dict(
            u0='sin(x[0])*cos(x[1])*cos(x[2])',
            u1='-cos(x[0])*sin(x[1])*cos(x[2])',
            u2='0',
            p='1./16.*(cos(2*x[0])+cos(2*x[1]))*(cos(2*x[2])+2)')))


def mesh(Nx, Ny, Nz, L, **params):
    mesh = BoxMesh(Point(-L / 2, -L / 2, -L / 2), Point(L / 2, L / 2, L / 2), Nx, Ny, Nz)
    print_mesh_information(mesh)
    return mesh


def near(x, y, tol=1e-12):
    return bool(abs(x - y) < tol)


def pre_boundary_condition(mesh, L, **NS_namespace):
    # Mark geometry
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    inlet = AutoSubDomain(lambda x, b: b)
    inlet.mark(boundary, 1)

    return dict(boundary=boundary)


class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        return bool((near(x[0], -pi) or near(x[1], -pi) or near(x[2], -pi)) and
                    (not (near(x[0], pi) or near(x[1], pi) or near(x[2], pi))) and on_boundary)

    def map(self, x, y):
        if near(x[0], pi) and near(x[1], pi) and near(x[2], pi):
            y[0] = x[0] - 2.0 * pi
            y[1] = x[1] - 2.0 * pi
            y[2] = x[2] - 2.0 * pi
        elif near(x[0], pi) and near(x[1], pi):
            y[0] = x[0] - 2.0 * pi
            y[1] = x[1] - 2.0 * pi
            y[2] = x[2]
        elif near(x[1], pi) and near(x[2], pi):
            y[0] = x[0]
            y[1] = x[1] - 2.0 * pi
            y[2] = x[2] - 2.0 * pi
        elif near(x[1], pi):
            y[0] = x[0]
            y[1] = x[1] - 2.0 * pi
            y[2] = x[2]
        elif near(x[0], pi) and near(x[2], pi):
            y[0] = x[0] - 2.0 * pi
            y[1] = x[1]
            y[2] = x[2] - 2.0 * pi
        elif near(x[0], pi):
            y[0] = x[0] - 2.0 * pi
            y[1] = x[1]
            y[2] = x[2]
        else:  # near(x[2], pi):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - 2.0 * pi


def pre_solve_hook(V, mesh, A0, T_G, L, t, newfolder, initial_fields_w, velocity_degree, u_components,
                   boundary, NS_expressions,
                   **NS_namespace):
    # Visualization files
    viz_p, viz_u = get_visualization_writers(newfolder, ["pressure", "velocity"])

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Mesh velocity conditions
    w0 = Expression((initial_fields_w["w0"]), element=V.ufl_element(), t=t, A=A0, T_G=T_G, L=L)
    w1 = Expression((initial_fields_w["w1"]), element=V.ufl_element(), t=t, A=A0, T_G=T_G, L=L)
    w2 = Expression((initial_fields_w["w2"]), element=V.ufl_element(), t=t, A=A0, T_G=T_G, L=L)

    NS_expressions["bc_w0"] = w0
    NS_expressions["bc_w1"] = w1
    NS_expressions["bc_w2"] = w2

    bc_mesh = dict((ui, []) for ui in u_components)
    bc0 = DirichletBC(V, w0, boundary, 1)
    bc1 = DirichletBC(V, w1, boundary, 1)
    bc2 = DirichletBC(V, w2, boundary, 1)

    bc_mesh["u0"] = [bc0]
    bc_mesh["u1"] = [bc1]
    bc_mesh["u2"] = [bc2]
    return dict(viz_p=viz_p, viz_u=viz_u, u_vec=u_vec, dof_map=dof_map, bc_mesh=bc_mesh, coordinates=coordinates)


def initialize(q_, q_1, q_2, VV, initial_fields, OasisFunction, **NS_namespace):
    for ui in q_:
        vv = OasisFunction(Expression((initial_fields[ui]),
                                      element=VV[ui].ufl_element()), VV[ui])
        vv()
        q_[ui].vector()[:] = vv.vector()[:]
        if not ui == 'p':
            q_1[ui].vector()[:] = q_[ui].vector()[:]
            q_2[ui].vector()[:] = q_[ui].vector()[:]


def update_boundary_conditions(t, NS_expressions, **NS_namespace):
    for key, value in NS_expressions.items():
        if "bc" in key:
            value.t = t


def temporal_hook(q_, save_solution_frequency, u_vec, viz_u, tstep, t, **NS_namespace):
    # Save velocity and pressure
    if tstep % save_solution_frequency == 0:
        assign(u_vec.sub(0), q_["u0"])
        assign(u_vec.sub(1), q_["u1"])
        assign(u_vec.sub(2), q_["u2"])

        viz_u.write(u_vec, t)
