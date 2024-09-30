import pickle
from os import getcwd

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers

comm = MPI.comm_world


def problem_parameters(
    commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace
):
    """
    Problem file for running CFD simulation for the moving Taylor-Green problem in 3D as described by Fehn et al.[1].
    The problem solves the N-S equations in the absence of body forces, and is commonly used to study transitional and
    turbulent flows. The problem initializes the solution at the two previous time steps, and applies periodic boundary
    condition on the domain walls in all coordinate directions. The domain boundaries are moving for the given mesh
    motion, while the fluid velocity is defined periodically in order to ensure consistency with the periodic boundary
    conditions. The main simulation parameters are the ampliture A, and period time duration T_G.

    [1] Fehn, N., Heinz, J., Wall, W. A., & Kronbichler, M. (2021). High-order arbitrary Lagrangian–Eulerian
    discontinuous Galerkin methods for the incompressible Navier–Stokes equations.
    Journal of Computational Physics, 430, 110040.
    """

    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(getcwd(), restart_folder)

        f = open(
            path.join(
                path.dirname(path.abspath(__file__)), restart_folder, "params.dat"
            ),
            "rb",
        )
        NS_parameters.update(pickle.load(f))
        NS_parameters["restart_folder"] = restart_folder
        globals().update(NS_parameters)

    else:
        # Override some problem specific parameters
        NS_parameters.update(
            # Problem specific parameters
            L=2 * np.pi,  # Mesh size
            A0=np.pi / 6,  # Amplitude
            T_G=20,  # Period time
            # Fluid parameters
            Re=1600,
            # Simulation parameters
            T=1,  # End time
            dt=0.005,  # Time step
            Nx=32,  # Resolution in the x-direction
            Ny=32,  # Resolution in the y-direction
            Nz=32,  # Resolution in the z-direction
            folder="results_moving_taylor_green_3d",
            # Oasis parameters
            max_iter=2,
            dynamic_mesh=True,
            save_solution_frequency=5e10,
            save_step=5,
            checkpoint=500,
            print_intermediate_info=100,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
        )

    NS_expressions.update(
        dict(
            constrained_domain=PeriodicDomain(),
            kin=np.zeros(1),
            initial_fields_w=dict(
                w0="2 * pi * A / T_G * cos(2 * pi * t / T_G)"
                + "* sin(2 * pi * (x[1] + L/2)/L) * sin(2*pi* (x[2] + L/2)/L)",
                w1="2*pi*A/T_G*cos(2*pi*t/T_G) * sin(2*pi * (x[0] + L/2)/L) * sin(2*pi* (x[2] + L/2)/L)",
                w2="2*pi*A/T_G*cos(2*pi*t/T_G) * sin(2*pi * (x[0] + L/2)/L) * sin(2*pi* (x[1] + L/2)/L)",
            ),
            initial_fields=dict(
                u0="sin(x[0])*cos(x[1])*cos(x[2])",
                u1="-cos(x[0])*sin(x[1])*cos(x[2])",
                u2="0",
                p="1./16.*(cos(2*x[0])+cos(2*x[1]))*(cos(2*x[2])+2)",
            ),
        )
    )


def mesh(Nx, Ny, Nz, L, dt, **params):
    mesh = BoxMesh(
        Point(-L / 2, -L / 2, -L / 2), Point(L / 2, L / 2, L / 2), Nx, Ny, Nz
    )

    print_mesh_information(mesh, dt, u_mean=1)
    return mesh


def near(x, y, tol=1e-12):
    return bool(abs(x - y) < tol)


def pre_boundary_condition(mesh, Re, **NS_namespace):
    # Mark geometry
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    inlet = AutoSubDomain(lambda x, b: b)
    inlet.mark(boundary, 1)

    # Update viscosity
    nu = 1 / Re
    return dict(boundary=boundary, nu=nu)


class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        return bool(
            (near(x[0], -pi) or near(x[1], -pi) or near(x[2], -pi))
            and (not (near(x[0], pi) or near(x[1], pi) or near(x[2], pi)))
            and on_boundary
        )

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


def pre_solve_hook(
    V,
    mesh,
    A0,
    T_G,
    L,
    t,
    newfolder,
    initial_fields_w,
    velocity_degree,
    u_components,
    boundary,
    NS_expressions,
    **NS_namespace,
):
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
    bc_mesh = dict((ui, []) for ui in u_components)
    for u, w in zip(["u0", "u1", "u2"], ["w0", "w1", "w2"]):
        wall_motion = Expression(
            (initial_fields_w[w]), element=V.ufl_element(), t=t, A=A0, T_G=T_G, L=L
        )
        NS_expressions[f"bc_{w}"] = wall_motion

        bc = DirichletBC(V, wall_motion, boundary, 1)
        bc_mesh[u] = [bc]

    return dict(
        viz_p=viz_p,
        viz_u=viz_u,
        u_vec=u_vec,
        dof_map=dof_map,
        bc_mesh=bc_mesh,
        coordinates=coordinates,
    )


def initialize(q_, q_1, q_2, VV, initial_fields, OasisFunction, **NS_namespace):
    for ui in q_:
        vv = OasisFunction(
            Expression((initial_fields[ui]), element=VV[ui].ufl_element()), VV[ui]
        )
        vv()
        q_[ui].vector()[:] = vv.vector()[:]
        if not ui == "p":
            q_1[ui].vector()[:] = q_[ui].vector()[:]
            q_2[ui].vector()[:] = q_[ui].vector()[:]


def update_boundary_conditions(t, NS_expressions, **NS_namespace):
    for key, value in NS_expressions.items():
        if "bc" in key:
            value.t = t


def temporal_hook(
    nu, L, dt, mesh, u_, save_solution_frequency, u_vec, viz_u, tstep, t, **NS_namespace
):
    # Save velocity and pressure
    if tstep % save_solution_frequency == 0:
        for i in range(3):
            assign(u_vec.sub(i), u_[i])

        viz_u.write(u_vec, t)

    # Compute mean velocity and Reynolds number at inlet
    if tstep % 10 == 0:
        h = mesh.hmin()
        compute_flow_quantities(
            u_,
            L,
            nu,
            mesh,
            t,
            tstep,
            dt,
            h,
            outlet_area=1,
            boundary=None,
            outlet_ids=[],
            inlet_ids=[],
            id_wall=0,
            period=1.0,
            newfolder=None,
            dynamic_mesh=False,
            write_to_file=False,
        )
