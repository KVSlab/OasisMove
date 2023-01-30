import json
import pickle
from pprint import pprint

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingAtriumCommon import *
from oasismove.problems.NSfracStep.MovingCommon import get_visualization_writers

set_log_level(50)


def problem_parameters(commandline_kwargs, NS_parameters, scalar_components, Schmidt, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
    else:
        backflow = bool(commandline_kwargs.get("backflow", True))
        backflow_beta = float(commandline_kwargs.get("backflow_beta", 0.5))
        # Override some problem specific parameters
        # Parameters are in mm and ms
        nu = 3.3018868e-3
        NS_parameters.update(
            # Atrium specific parameters
            flow_rate_type=None,  # Type of flow rate profile - [AF,SR,Rigid]
            profile_power=2,  # Degree of velocity profile
            movement_at_inlet_and_outlet=False,  # Adds movement to inlet and outlet boundaries
            dynamic_mesh=True,  # Run moving mesh simulation
            compute_velocity_and_pressure=True,  # Only solve mesh equations
            # Fluid parameters
            # Viscosity [nu: 0.0035 Pa-s / 1060 kg/m^3 = 3.3018868E-6 m^2/s == 3.3018868E-3 mm^2/ms]
            nu=nu,
            # Geometry parameters
            id_in=[],  # Inlet boundary ID
            id_out=[],  # Outlet boundary IDs
            # Simulation parameters
            cycle=1000,  # Run simulation for 1 cardiac cycles [ms]
            no_of_cycles=2,  # Number of cycles
            dt=1,  # dt=0.01 <=> 10 000 steps per cycle [ms]
            checkpoint=20,  # Overwrite solution in Checkpoint folder each checkpoint
            save_step=1e10,  # Store solution each save_step
            # Frequencies to save data
            save_solution_tstep=5,  # Store velocity at this freq
            save_flow_metrics_tstep=5,  # Store flow metrics at this freq
            save_volume_tstep=5e5,  # Compute and store volume at this freq
            print_intermediate_info=10,
            folder="results_moving_atrium",
            mesh_path=commandline_kwargs["mesh_path"],
            # Solver parameters
            max_iter=2,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True,
            krylov_solvers=dict(monitor_convergence=False)
        )

        if backflow:
            NS_parameters.update(
                # Backflow
                backflow_facets=[],  # Backflow stabilization at outlet
                backflow_beta=backflow_beta
            )

    mesh_file = NS_parameters["mesh_path"].split("/")[-1]
    case_name = mesh_file.split(".")[0]
    NS_parameters["folder"] = path.join(NS_parameters["folder"], case_name)

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Starting simulation for MovingAtrium.py ===")
        print("Running with the following parameters:")
        pprint(NS_parameters)


def scalar_source(**NS_namespace):
    return dict(blood=Constant(1.0))


def mesh(mesh_path, **NS_namespace):
    # Read mesh and print mesh information
    atrium_mesh = Mesh(mesh_path)
    print_mesh_information(atrium_mesh)

    return atrium_mesh


def pre_boundary_condition(mesh, mesh_path, id_out, id_in, no_of_cycles, cycle, flow_rate_type, backflow_facets,
                           NS_parameters, **NS_namespace):
    # Variables needed during the simulation
    D = mesh.geometry().dim()
    boundary = MeshFunction("size_t", mesh, D - 1, mesh.domains())
    ds_new = Measure("ds", domain=mesh, subdomain_data=boundary)

    # Get IDs for inlet(s) and outlet(s) and set backflow at outlet
    info_path = mesh_path.split(".")[0] + "_info.json"
    with open(info_path) as f:
        info = json.load(f)

    id_wall = 0
    id_in[:] = info['inlet_ids']
    id_out[:] = info['outlet_id']
    # FIXME: Not working in parallel
    # backflow_facets[:] = info['outlet_id']

    # Find corresponding areas
    outlet_area = info['outlet_area']
    area_total = 0
    for ID in id_in:
        area_total += assemble(Constant(1.0) * ds_new(ID))

    D_mitral = np.sqrt(4 * outlet_area / np.pi)

    T = no_of_cycles * cycle

    # Load flow rate
    s_pv = 1E2
    if flow_rate_type == "AFMoving":
        flow_rate_data = np.loadtxt(path.join(mesh_path.rsplit("/", 1)[0], "../Q_AF_moving.txt"))
        weight = np.ones(len(flow_rate_data))
        weight[-5:] = 0.01
    elif flow_rate_type == "AFRigid":
        flow_rate_data = np.loadtxt(path.join(mesh_path.rsplit("/", 1)[0], "../Q_AF_rigid.txt"))
        weight = np.ones(len(flow_rate_data))
        weight[-20:] = 0.01
    elif flow_rate_type == "SRMoving":
        flow_rate_data = np.loadtxt(path.join(mesh_path.rsplit("/", 1)[0], "../Q_SR_moving.txt"))
        weight = np.ones(len(flow_rate_data))
    elif flow_rate_type == "SRRigid":
        flow_rate_data = np.loadtxt(path.join(mesh_path.rsplit("/", 1)[0], "../Q_SR_rigid.txt"))
        weight = np.ones(len(flow_rate_data))
        weight[170:180] = 0.0001
    elif flow_rate_type == "LA5":
        flow_rate_data = np.loadtxt(path.join(mesh_path.rsplit("/", 1)[0], "LA5_flow_rate.txt"))
        weight = np.ones(len(flow_rate_data))
        s_pv = 50
    else:
        print("--- No flow rate type supplied. Exiting..")
        exit()

    # Spline flow rate data
    t_e = np.linspace(0, cycle, len(flow_rate_data))
    flow_rate = splrep(t_e, flow_rate_data, s=s_pv, per=True, w=weight)

    return dict(boundary=boundary, ds_new=ds_new, area_total=area_total, flow_rate=flow_rate, id_wall=id_wall, T=T,
                outlet_area=outlet_area, D_mitral=D_mitral)


def create_bcs(NS_expressions, cycle, x_, boundary, V, Q, mesh, mesh_path, id_in, id_out, id_wall, tstep, dt,
               dynamic_mesh, area_total, flow_rate, t, profile_power, flow_rate_type, **NS_namespace):
    # Create inlet boundary conditions
    coords = ['x', 'y', 'z']
    N = int(cycle / dt)
    bc_inlets = {}
    for i, ID in enumerate(id_in):
        tmp_area, tmp_center, tmp_radius, tmp_normal = compute_boundary_geometry_acrn(mesh, id_in[i], boundary)
        inlet = []
        for normal_component in tmp_normal:
            R2 = tmp_area / np.pi  # R**2

            inlet_expression = InletParabolic(tstep, dt, N, normal_component, flow_rate, tmp_center, R2, tmp_area,
                                              area_total, cycle, profile_power, element=V.ufl_element())
            inlet.append(inlet_expression)

        NS_expressions[ID] = inlet
        bc_inlet = [DirichletBC(V, inlet[i], boundary, ID) for i in range(len(coords))]
        bc_inlets[ID] = bc_inlet

    # Set outlet boundary conditions
    bc_p = []
    for i, ID in enumerate(id_out):
        bc = DirichletBC(Q, Constant(0), boundary, ID)
        bc_p.append(bc)
        NS_expressions['P'] = bc

    # Set wall motion if not rigid domain, i.e. moving boundaries
    if dynamic_mesh:
        # Get mesh motion
        print("Loading displacement points")
        points = np.load(mesh_path.split(".")[0] + "_points.np", allow_pickle=True)

        # Define wall movement
        wall_counter = Surface_counter(points, cycle, element=V.ufl_element())
        bc_tmp = DirichletBC(V, wall_counter, boundary, id_wall)
        bc_tmp.apply(x_["u0"])
        x_["u0"].zero()

        wall_counter_max = wall_counter.counter
        wall_motion = wall_counter.get_motion()

        # Define outlet movement
        outlet_counter = Surface_counter(points, cycle, element=V.ufl_element())
        bc_tmp = DirichletBC(V, outlet_counter, boundary, id_out[0])
        bc_tmp.apply(x_["u0"])
        x_["u0"].zero()

        outlet_counter_max = outlet_counter.counter
        outlet_motion = outlet_counter.get_motion()

        # Remove explicitly from memory.
        del wall_counter
        del outlet_counter

        for i, coord in enumerate(coords):
            wall_ = Wall_motion(t, wall_motion, wall_counter_max, i, cycle, element=V.ufl_element())
            outlet_ = Wall_motion(t, outlet_motion, outlet_counter_max, i, cycle, element=V.ufl_element())
            NS_expressions["wall_%s" % coord] = wall_
            NS_expressions["outlet_%s" % coord] = outlet_

        for i, ID in enumerate(id_in):
            # Define inlet movement
            inlet_counter = Surface_counter(points, cycle, element=V.ufl_element())
            bc_tmp = DirichletBC(V, inlet_counter, boundary, ID)
            bc_tmp.apply(x_["u0"])
            x_["u0"].zero()

            inlet_counter_max = inlet_counter.counter
            inlet_motion = inlet_counter.get_motion()

            for j, coord in enumerate(coords):
                inlet_ = Wall_motion(t, inlet_motion, inlet_counter_max, j, cycle, element=V.ufl_element())
                NS_expressions["inlet_%i_%s" % (ID, coord)] = inlet_

            # Remove explicitly from memory.
            del inlet_counter

        del points

    #  Fluid velocity at walls
    if dynamic_mesh:
        bcu_wall = [DirichletBC(V, NS_expressions["wall_%s" % coord], boundary, id_wall) for coord in coords]
    else:
        bcu_wall = [DirichletBC(V, Constant(0.0), boundary, id_wall)] * len(coords)

    # Create lists with all boundary conditions
    bc_u0 = []
    bc_u1 = []
    bc_u2 = []
    bc_u0.append(bcu_wall[0])
    bc_u1.append(bcu_wall[1])
    bc_u2.append(bcu_wall[2])
    for ID in id_in:
        bc_u0.append(bc_inlets[ID][0])
        bc_u1.append(bc_inlets[ID][1])
        bc_u2.append(bc_inlets[ID][2])

    return dict(u0=bc_u0, u1=bc_u1, u2=bc_u2, p=bc_p)


def pre_solve_hook(mesh, V, id_in, id_out, id_wall, newfolder, u_components, velocity_degree, x_, boundary,
                   dynamic_mesh, restart_folder, movement_at_inlet_and_outlet, **NS_namespace):
    # Create point for evaluation
    n = FacetNormal(mesh)
    quantities = ["pressure", "velocity"]
    viz_p, viz_u = get_visualization_writers(newfolder, quantities)

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    # Extract dof map and coordinates
    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    print("Setup mesh velocity")

    if dynamic_mesh:
        # Set wall motion BCS
        bc_mesh = dict((ui, []) for ui in u_components)

        bcu_in_x = []
        bcu_in_y = []
        bcu_in_z = []
        # Add mesh movement to inlet and outlet
        if movement_at_inlet_and_outlet:
            bc_out_x = DirichletBC(V, NS_expressions["outlet_x"], boundary, id_out[0])
            bc_out_y = DirichletBC(V, NS_expressions["outlet_y"], boundary, id_out[0])
            bc_out_z = DirichletBC(V, NS_expressions["outlet_z"], boundary, id_out[0])
            for i, ID in enumerate(id_in):
                bcu_x = DirichletBC(V, NS_expressions["inlet_%i_x" % ID], boundary, ID)
                bcu_y = DirichletBC(V, NS_expressions["inlet_%i_y" % ID], boundary, ID)
                bcu_z = DirichletBC(V, NS_expressions["inlet_%i_z" % ID], boundary, ID)
                bcu_in_x.append(bcu_x)
                bcu_in_y.append(bcu_y)
                bcu_in_z.append(bcu_z)
        else:
            noslip = Constant(0.0)

            bc_out_x = DirichletBC(V, noslip, boundary, id_out[0])
            bc_out_y = DirichletBC(V, noslip, boundary, id_out[0])
            bc_out_z = DirichletBC(V, noslip, boundary, id_out[0])
            for i, ID in enumerate(id_in):
                bcu = DirichletBC(V, noslip, boundary, ID)
                bcu_in_x.append(bcu)
                bcu_in_y.append(bcu)
                bcu_in_z.append(bcu)

        # Add wall movement to wall
        bcu_wall_x = DirichletBC(V, NS_expressions["wall_x"], boundary, id_wall)
        bcu_wall_y = DirichletBC(V, NS_expressions["wall_y"], boundary, id_wall)
        bcu_wall_z = DirichletBC(V, NS_expressions["wall_z"], boundary, id_wall)

        bc_mesh["u0"] = [bcu_wall_x] + [bc_out_x] + bcu_in_x
        bc_mesh["u1"] = [bcu_wall_y] + [bc_out_y] + bcu_in_y
        bc_mesh["u2"] = [bcu_wall_z] + [bc_out_z] + bcu_in_z
    else:
        bc_mesh = {}

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(newfolder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    # Save mesh as HDF5 file
    with HDF5File(MPI.comm_world, files["mesh"], "w") as mesh_file:
        mesh_file.write(mesh, "mesh")

    return dict(dof_map=dof_map, coordinates=coordinates, boundary=boundary, viz_u=viz_u, viz_p=viz_p, u_vec=u_vec, n=n,
                bc_mesh=bc_mesh)


def update_boundary_conditions(t, dynamic_mesh, NS_expressions, id_in, tstep, **NS_namespace):
    # Update inlet condition
    for ID in id_in:
        for i in [0, 1, 2]:
            NS_expressions[ID][i].update(tstep)

    if dynamic_mesh:
        # Update wall motion BCs
        for coord in ["x", "y", "z"]:
            NS_expressions["wall_{}".format(coord)].t = t
            NS_expressions["outlet_{}".format(coord)].t = t
            for i, ID in enumerate(id_in):
                NS_expressions["inlet_%i_%s" % (ID, coord)].t = t


def velocity_tentative_hook(mesh, boundary, u_ab, x_1, b, A, ui, u, v, backflow_facets, backflow_beta, **NS_namespace):
    add_backflow_stabilization(A, b, backflow_beta, backflow_facets, boundary, mesh, u, u_ab, ui, v, x_1)


def temporal_hook(u_, mesh, newfolder, boundary, tstep, t, nu, dt, u_vec, viz_u, id_in, id_out, id_wall, outlet_area,
                  cycle, D_mitral, p_, viz_p, save_solution_tstep, save_flow_metrics_tstep, save_volume_tstep,
                  **NS_namespace):
    if tstep % save_volume_tstep == 0:
        compute_volume(mesh, t, newfolder)

    if tstep % save_solution_tstep == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])
        assign(u_vec.sub(2), u_[2])

        viz_u.write(u_vec, t)

    if tstep % save_flow_metrics_tstep == 0:
        h = mesh.hmin()
        compute_flow_quantities(u_, D_mitral, nu, mesh, t, tstep, dt, h, outlet_area, boundary, id_out, id_in, id_wall,
                                period=cycle, newfolder=newfolder, dynamic_mesh=False, write_to_file=True)
