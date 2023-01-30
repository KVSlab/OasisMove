import os
import pickle

from IPython import embed
from scipy.interpolate import splrep, splev

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingAtriumCommon import Surface_counter
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers


def problem_parameters(commandline_kwargs, NS_parameters, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(os.getcwd(), restart_folder)
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        # Override some problem specific parameters
        Re = 3475.59
        d_mitral = 0.42857
        u_mean = 15.516
        nu = 0.00191
        cardiac_cycle = 0.9
        number_of_cycles = 1
        Re_computed = u_mean * d_mitral / nu
        if MPI.rank(MPI.comm_world) == 0:
            print("Real RE={} | Computed RE={}".format(Re, Re_computed))
        NS_parameters.update(
            # Enable backflow stabilization
            backflow_facets=[2],
            backflow_beta=1.0,
            # Mesh parameters
            id_wall=0,
            id_aorta=1,
            id_mitral=2,
            d_mitral=d_mitral,
            # Simulation parameters
            u_mean=u_mean,
            number_of_cycles=number_of_cycles,
            cardiac_cycle=cardiac_cycle,  # s
            T=cardiac_cycle * number_of_cycles,
            dt=0.9/500,
            nu=nu,  # cm^2/s
            dynamic_mesh=True,
            checkpoint=500,
            save_step_problem=1,
            save_step_problem_h5=20E10,
            save_flow_metrics_tstep=40E10,
            folder="results_moving_ventricle",
            mesh_path="src/oasismove/mesh/MovingVentricle3D/Ventricle3D/wall_mesh_003.xml.gz",
            velocity_degree=1,
            pressure_degree=1,
            print_intermediate_info=20,
            max_iter=2,
            use_krylov_solvers=True)


def mesh(mesh_path, u_mean, dt, velocity_degree, **NS_namespace):
    # Read mesh and print mesh information
    mesh = Mesh(mesh_path)
    print_mesh_information(mesh)

    # Compute mesh information and CFL
    V = VectorFunctionSpace(mesh, "CG", velocity_degree)
    dofs = V.dim()
    dx = mesh.hmin()
    if MPI.rank(MPI.comm_world) == 0:
        print("N_cells/CPU = {} | CFL~{} | dx={} | dt = {} | dofs={}".format(mesh.num_cells(), u_mean * dt / dx, dx, dt,
                                                                         dofs))
    return mesh


def pre_boundary_condition(mesh, id_aorta, id_mitral, **NS_namespace):
    # Variables needed during the simulation
    D = mesh.geometry().dim()
    boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
    ds = Measure("ds", domain=mesh, subdomain_data=boundary)
    ds_out = ds(id_aorta)
    ds_in = ds(id_mitral)
    area_outlet = assemble(Constant(1.0) * ds_out)
    area_inlet = assemble(Constant(1.0) * ds_in)

    hmin = mesh.hmin()
    return dict(boundary=boundary, D=D, area_outlet=area_outlet, area_inlet=area_inlet, hmin=hmin)


class Wall_motion(UserExpression):
    def __init__(self, t, motion, max_counter, direction, cycle, **kwargs):
        self.t = t
        self.motion = motion
        self.max_counter = max_counter
        self.counter = -1
        self.direction = direction
        self.cycle = cycle
        super().__init__(**kwargs)

    def eval(self, values, x):
        self.counter += 1
        # No motion for inlet/outlet
        values[:] = splev(self.t % self.cycle, self.motion[self.counter][self.direction], der=1)
        if self.counter == self.max_counter:
            self.counter = -1


class InletPlug(UserExpression):
    def __init__(self, tstep, dt, N, f, L, cycle, **kwargs):
        self.tstep = tstep
        self.dt = dt
        self.N = N
        self.L = L
        self.f = f
        self.U0 = 0
        self.cycle = cycle
        super().__init__(**kwargs)

    def update(self, tstep):
        self.tstep = tstep
        tstep = self.tstep % self.N
        self.U0 = - splev(tstep * self.dt % self.cycle, self.f) / self.L

    def eval(self, value, _):
        value[:] = self.U0


def create_bcs(mesh_path, V, Q, NS_expressions, boundary, area_inlet, id_wall, id_aorta, id_mitral, t, x_,
               cardiac_cycle, dt, tstep, **NS_namespace):
    N = int(cardiac_cycle / dt)
    # Load and spline flow rate
    flow_rate_path = mesh_path.split(".xml")[0] + "_flowrate_moving.txt"
    s_m = 0.6E-2
    t_values, flow_rate = np.loadtxt(flow_rate_path).T
    inlet_profile = splrep(t_values, flow_rate, s=s_m, per=True)
    inlet_mitral = InletPlug(tstep, dt, N, inlet_profile, area_inlet, cardiac_cycle, element=V.ufl_element())
    NS_expressions["inlet_mitral"] = inlet_mitral

    # Set inlet flow profile at mitral valve
    bc_u_x = DirichletBC(V, 0, boundary, id_mitral)
    bc_u_y = DirichletBC(V, 0, boundary, id_mitral)
    bc_u_z = DirichletBC(V, inlet_mitral, boundary, id_mitral)

    # Pressure at aortic valve
    bc_p = DirichletBC(Q, Constant(0.0), boundary, id_aorta)

    # Read points for wall surface
    # Load displaced points
    if MPI.rank(MPI.comm_world) == 0:
        print("Loading displacement points")
    points = np.load(mesh_path.split(".")[0] + ".np", allow_pickle=True)

    # Get mesh motion
    wall_counter = Surface_counter(points, cardiac_cycle, element=V.ufl_element())
    bc_tmp = DirichletBC(V, wall_counter, boundary, id_wall)
    bc_tmp.apply(x_["u0"])
    x_["u0"].zero()

    counter_max = wall_counter.counter
    motion = wall_counter.get_motion()

    # Remove points explicitly from memory.
    del wall_counter
    del points

    # Expressions
    wall0 = Wall_motion(t, motion, counter_max, 0, cardiac_cycle, element=V.ufl_element())
    wall1 = Wall_motion(t, motion, counter_max, 1, cardiac_cycle, element=V.ufl_element())
    wall2 = Wall_motion(t, motion, counter_max, 2, cardiac_cycle, element=V.ufl_element())

    # Store for later
    NS_expressions["expressions_wall"] = [wall0, wall1, wall2]

    # BCs wall
    wall_bc_0 = DirichletBC(V, wall0, boundary, id_wall)
    wall_bc_1 = DirichletBC(V, wall1, boundary, id_wall)
    wall_bc_2 = DirichletBC(V, wall2, boundary, id_wall)
    return dict(p=[bc_p],
                u0=[wall_bc_0, bc_u_x],
                u1=[wall_bc_1, bc_u_y],
                u2=[wall_bc_2, bc_u_z])


def pre_solve_hook(mesh, newfolder, V, velocity_degree, u_components, NS_expressions, boundary,
                   id_mitral,
                   id_aorta, id_wall, restart_folder, **NS_namespace):
    # Visualization files
    viz_files = get_visualization_writers(newfolder, ['velocity'])

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    # Extract dof map and coordinates
    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Set wall motion BCs
    bc_mesh = dict((ui, []) for ui in u_components)
    bc_in = [DirichletBC(V, Constant(0.0), boundary, id_mitral)]
    bc_out = [DirichletBC(V, Constant(0.0), boundary, id_aorta)]

    [wall0, wall1, wall2] = NS_expressions["expressions_wall"]
    wall_bc_0 = DirichletBC(V, wall0, boundary, id_wall)
    wall_bc_1 = DirichletBC(V, wall1, boundary, id_wall)
    wall_bc_2 = DirichletBC(V, wall2, boundary, id_wall)

    bc_mesh["u0"] = [wall_bc_0] + bc_out + bc_in
    bc_mesh["u1"] = [wall_bc_1] + bc_out + bc_in
    bc_mesh["u2"] = [wall_bc_2] + bc_out + bc_in

    # For writing solution
    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(newfolder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    # Save mesh as HDF5 file for post processing
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")

    # Create vector function for storing velocity
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_mean = Function(Vv)
    u_mean0 = Function(V)
    u_mean1 = Function(V)
    u_mean2 = Function(V)
    U = Function(Vv)

    # Tstep when solutions for post processing should start being saved
    save_solution_at_tstep = 0

    return dict(viz_files=viz_files, coordinates=coordinates,
                dof_map=dof_map, bc_mesh=bc_mesh, u_vec=u_vec, u_mean=u_mean, u_mean0=u_mean0,
                save_solution_at_tstep=save_solution_at_tstep, u_mean1=u_mean1,
                u_mean2=u_mean2, U=U)


# def velocity_tentative_hook(mesh, boundary, u_ab, x_1, b, A, ui, u, v, backflow_facets, backflow_beta, **NS_namespace):
#     add_backflow_stabilization(A, b, backflow_beta, backflow_facets, boundary, mesh, u, u_ab, ui, v, x_1)


def update_prescribed_motion(t, u_components, NS_expressions, tstep, **NS_namespace):
    # Update wall motion BCs
    for i, ui in enumerate(u_components):
        NS_expressions["expressions_wall"][i].t = t

    for key, value in NS_expressions.items():
        if "inlet" in key:
            value.update(tstep)


def temporal_hook(d_mitral, save_flow_metrics_tstep, id_mitral, id_aorta, nu, boundary, newfolder, mesh, dt, U,
                  save_step_problem, u_vec, viz_files, t, tstep, u_, id_wall, save_step_problem_h5, NS_parameters,
                  cardiac_cycle, u_mean0, u_mean1, u_mean2, area_outlet, **NS_namespace):
    if tstep % save_step_problem == 0:
        print("Saving")
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])
        assign(u_vec.sub(2), u_[2])

        viz_files[0].write(u_vec, t)
        print("Done")

    if tstep % save_step_problem_h5 == 0:
        assign(U.sub(0), u_[0])
        assign(U.sub(1), u_[1])
        assign(U.sub(2), u_[2])

        # Get save paths
        files = NS_parameters['files']
        u_path = files['u']
        file_mode = "w" if not path.exists(u_path) else "a"

        # Save velocity
        viz_u = HDF5File(MPI.comm_world, u_path, file_mode=file_mode)
        viz_u.write(U, "/velocity", tstep)
        viz_u.close()

        # Accumulate velocity
        u_mean0.vector().axpy(1, u_[0].vector())
        u_mean1.vector().axpy(1, u_[1].vector())
        u_mean2.vector().axpy(1, u_[2].vector())

    if tstep % save_flow_metrics_tstep == 0:
        h = mesh.hmin()
        compute_flow_quantities(u_, d_mitral, nu, mesh, t, tstep, dt, h, area_outlet, boundary, [id_aorta], [id_mitral],
                                id_wall, period=cardiac_cycle, newfolder=newfolder, dynamic_mesh=False,
                                write_to_file=True)


def get_file_paths(folder):
    # Create folder where data and solutions (velocity, mesh, pressure) is stored
    common_path = path.join(folder, "Solutions")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path):
            os.makedirs(common_path)

    file_p = path.join(common_path, "p.h5")
    file_u = path.join(common_path, "u.h5")
    file_u_mean = path.join(common_path, "u_mean.h5")
    file_mesh = path.join(common_path, "mesh.h5")
    files = {"u": file_u, "u_mean": file_u_mean, "p": file_p, "mesh": file_mesh}

    return files
