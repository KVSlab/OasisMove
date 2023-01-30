import os
import pickle

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_files


def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    """
    Problem file for running CFD simulation for the oscillating cylinder in a rectangular 2D domain, as described
    by Blackburn and Henderson [1].The cylinder is prescribed an oscillatory motion and is placed in a free stream,
    with a diameter of D cm. The kinetmatic viscosity is computed from the free stream velocity of 1 m/s for a Reynolds
    number of 500, which is well above the critical value for vortex shedding. Moreover, the oscillation is mainly
    controlled by the amplitude ratio A_ratio, the Strouhal number St, and the frequency ratio F.

    [1] Blackburn, H. M., & Henderson, R. D. (1999). A study of two-dimensional flow past an oscillating cylinder.
    Journal of Fluid Mechanics, 385, 255-286.
    """

    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(os.getcwd(), restart_folder)

        f = open(path.join(path.dirname(path.abspath(__file__)), restart_folder,
                           'params.dat'), 'r')
        NS_parameters.update(pickle.load(f))
        NS_parameters['T'] = 10
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)

    else:
        T = 5  # End time
        Re = 500  # Reynolds number
        D = 0.1  # Cylinder diameter
        U = 1.0  # Free stream velocity
        nu = U * D / Re
        NS_parameters.update(
            # Geometrical parameters
            diam=D,  # Diameter
            nu=nu,  # Kinematic viscosity
            U=U,  # Free-stream flow speed
            A_ratio=0.25,  # Amplitude ratio
            St=0.2280,  # Strouhal number
            F=1.0,  # Frequency ratio
            # Simulation parameters
            T=T,
            dt=0.000125,
            checkpoint=500,
            save_solution_frequency=5,
            print_intermediate_info=100,
            velocity_degree=1,
            pressure_degree=1,
            mesh_path=commandline_kwargs["mesh_path"],
            folder="results_moving_cylinder",
            use_krylov_solvers=True,
            max_iter=2,
            max_error=1E-8)


def mesh(mesh_path, dt, U, **NS_namespace):
    # Import mesh
    mesh = Mesh()
    with XDMFFile(MPI.comm_world, mesh_path) as infile:
        infile.read(mesh)

    print("Mesh info: N_cells = {} |  dx={} | dt = {}".format(mesh.num_cells(), mesh.hmin(), dt))
    print("CFL = {} | U = {}".format(U * dt / mesh.hmin(), U))
    return mesh


def pre_boundary_condition(diam, mesh, **NS_namespace):
    # Reference domain from Blackburn & Henderson (1999)
    # Cylinder centered in (0,0)
    H = 30 * diam / 2  # Height
    L1 = -10 * diam  # Length
    L2 = 52 * diam  # Length

    # Mark geometry
    inlet = AutoSubDomain(lambda x, b: b and x[0] <= L1 + DOLFIN_EPS)
    walls = AutoSubDomain(lambda x, b: b and (near(x[1], -H) or near(x[1], H)))
    circle = AutoSubDomain(lambda x, b: b and (-H / 2 <= x[1] <= H / 2) and (L1 / 2 <= x[0] <= L2 / 2))
    outlet = AutoSubDomain(lambda x, b: b and (x[0] > L2 - DOLFIN_EPS * 1000))

    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    inlet.mark(boundary, 1)
    walls.mark(boundary, 2)
    circle.mark(boundary, 3)
    outlet.mark(boundary, 4)

    return dict(boundary=boundary)


def create_bcs(V, Q, U, St, F, diam, A_ratio, sys_comp, boundary, NS_expressions, **NS_namespace):
    info_red("Creating boundary conditions")

    f_v = St * U / diam  # Fixed-cylinder vortex shredding frequency
    f_o = F * f_v  # Frequency of harmonic oscillation
    y_max = A_ratio * diam  # Max displacement (Amplitude)
    print("Frequency is %.4f" % f_o)
    print("Amplitude is %.4f " % y_max)

    NS_expressions["circle_x"] = Constant(0)
    NS_expressions["circle_y"] = Expression('2 * pi * f_o * y_max* cos(2 * pi * f_o * t)', degree=2, t=0, y_max=y_max,
                                            f_o=f_o)

    bcu_in_x = DirichletBC(V, Constant(U), boundary, 1)
    bcu_in_y = DirichletBC(V, Constant(0), boundary, 1)

    bcu_wall_x = DirichletBC(V, Constant(U), boundary, 2)
    bcu_wall_y = DirichletBC(V, Constant(0), boundary, 2)

    bcu_circle_x = DirichletBC(V, NS_expressions["circle_x"], boundary, 3)
    bcu_circle_y = DirichletBC(V, NS_expressions["circle_y"], boundary, 3)

    bcp_out = DirichletBC(Q, Constant(0), boundary, 4)

    bcs = dict((ui, []) for ui in sys_comp)
    bcs['u0'] = [bcu_circle_x, bcu_in_x, bcu_wall_x]
    bcs['u1'] = [bcu_circle_y, bcu_in_y, bcu_wall_y]
    bcs["p"] = [bcp_out]

    return bcs


def pre_solve_hook(V, p_, u_, velocity_degree, nu, mesh, newfolder, u_components, boundary, **NS_namespace):
    # Visualization files
    viz_p, viz_u = get_visualization_files(newfolder)

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Inlet, walls
    rigid_bc_in = DirichletBC(V, Constant(0), boundary, 1)
    rigid_bc_walls = DirichletBC(V, Constant(0), boundary, 2)
    circle_bc_x = DirichletBC(V, NS_expressions["circle_x"], boundary, 3)
    circle_bc_y = DirichletBC(V, NS_expressions["circle_y"], boundary, 3)
    rigid_bc_out = DirichletBC(V, Constant(0), boundary, 4)

    bc_mesh = dict((ui, []) for ui in u_components)
    rigid_bc = [rigid_bc_in, rigid_bc_walls, rigid_bc_out]
    bc_mesh["u0"] = [circle_bc_x] + rigid_bc
    bc_mesh["u1"] = [circle_bc_y] + rigid_bc

    ds = Measure("ds", subdomain_data=boundary)

    R = VectorFunctionSpace(mesh, 'R', 0)
    c = TestFunction(R)
    n = FacetNormal(mesh)
    tau = -p_ * Identity(2) + nu * (grad(u_) + grad(u_).T)
    forces = dot(dot(tau, n), c) * ds(3)

    return dict(dof_map=dof_map, viz_p=viz_p, viz_u=viz_u, u_vec=u_vec, forces=forces, coordinates=coordinates)


def update_boundary_conditions(t, NS_expressions, **NS_namespace):
    # Update time
    NS_expressions["circle_y"].t = t


def temporal_hook(t, St, F, A_ratio, tstep, save_solution_frequency, forces, q_, U, diam, viz_u, viz_p, u_vec,
                  newfolder,
                  **NS_namespace):
    # Save fluid velocity and pressure solution
    if tstep % save_solution_frequency == 0:
        assign(u_vec.sub(0), q_["u0"])
        assign(u_vec.sub(1), q_["u1"])

        viz_u.write(u_vec, t)
        viz_p.write(q_["p"], t)

    # Compute drag and lift coefficients
    rho = 1000
    factor = 0.5 * rho * U ** 2 * diam
    Dr = assemble(forces).get_local() * rho  # Times constant fluid density

    f_v = St * U / diam  # Fixed-cylinder vortex shredding frequency
    f_o = F * f_v  # Frequency of harmonic oscillation
    y_max = A_ratio * diam  # Max displacement (Amplitude)

    # Compute pressure values
    rho = 1000
    dy_current = y_max * np.sin(2 * np.pi * f_o * t)
    p_0 = p_180 = 0
    try:
        p_0 = q_["p"](-diam / 2, dy_current)
    except KeyError:
        pass
    try:
        p_180 = q_["p"](diam / 2, dy_current)
    except KeyError:
        pass

    C_pb = 1 + 2 * (p_180 - p_0) / (rho * U ** 2)

    # Write forces to file
    fluid_properties = np.asarray([round(t, 4), Dr[0], Dr[0] / factor, Dr[1], Dr[1] / factor, C_pb])
    forces_path = path.join(newfolder, "Solutions", "forces.txt")
    if not path.isdir(path.join(newfolder, "Solutions")):
        os.mkdir(path.join(newfolder, "Solutions"))
    with open(forces_path, 'a') as filename:
        filename.write("{:.4f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} \n".format(*fluid_properties))
