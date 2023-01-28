from pprint import pprint

from .MovingCommon import get_visualization_writers, get_coordinate_map
from ..NSfracStep import *


# Override some problem specific parameters
def problem_parameters(NS_parameters, NS_expressions, commandline_kwargs, **NS_namespace):
    L = 1.0
    T = 1.0
    T_G = 4 * T
    T = 4
    dt = 5 * 10 ** (-2)
    A = 0.08  # Amplitude
    try:
        unstructured = commandline_kwargs["umesh"]
    except KeyError:
        unstructured = False

    NS_parameters.update(
        dynamic_mesh=True,
        unstructured=unstructured,
        nu=0.025,
        T=T,
        A0=A,
        T_G=T_G,
        L=L,
        dt=dt,
        Nx=20, Ny=20,
        folder="results_moving_vortex",
        plot_interval=1000,
        save_step=10000,
        checkpoint=10000,
        print_intermediate_info=10000,
        compute_error=1,
        use_krylov_solvers=True,
        velocity_degree=1,
        pressure_degree=1,
        max_iter=2,
        krylov_report=False)

    NS_parameters['krylov_solvers'] = {'monitor_convergence': False,
                                       'report': False,
                                       'relative_tolerance': 1e-8,
                                       'absolute_tolerance': 1e-8}

    NS_expressions.update(dict(
        exact_fields=dict(
            u0='-sin(2*pi*x[1])*exp(-4*pi*pi*nu*t)',
            u1=' sin(2*pi*x[0])*exp(-4*pi*pi*nu*t)',
            p=' -cos(2*pi*x[0])*cos(2*pi*x[1])*exp(-8*pi*pi*nu*t)'),
        initial_fields=dict(
            u0='-sin(2*pi*x[1])*exp(-4*pi*pi*nu*t)',
            u1=' sin(2*pi*x[0])*exp(-4*pi*pi*nu*t)',
            p=' -cos(2*pi*x[0])*cos(2*pi*x[1])*exp(-8*pi*pi*nu*t)'),
        initial_fields_w=dict(
            w0='A*2*pi / T_G * cos(2*pi*t/T_G)*sin(2*pi*(x[1] + L/2)/L)',
            w1='A*2*pi / T_G * cos(2*pi*t/T_G)*sin(2*pi*(x[0] + L/2)/L)'),
        total_error=np.zeros(3)))

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Starting simulation for MovingVortex.py ===")
        print("Running with the following parameters:")
        pprint(NS_parameters)


def mesh(Nx, Ny, L, dt, **params):
    mesh = RectangleMesh(Point(-L / 2, -L / 2), Point(L / 2, L / 2), Nx, Ny)

    print("Mesh info: N_cells = {} |  dx={} | dt = {}".format(mesh.num_cells(), mesh.hmin(), dt))
    return mesh


def pre_boundary_condition(mesh, **NS_namespace):
    # Mark geometry
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    inlet = AutoSubDomain(lambda x, b: b)

    inlet.mark(boundary, 1)
    return dict(boundary=boundary)


def create_bcs(V, Q, t, dt, nu, sys_comp, boundary, initial_fields, **NS_namespace):
    for i, ui in enumerate(sys_comp):
        if 'IPCS' in NS_parameters['solver']:
            deltat = dt / 2. if ui == 'p' else 0.
        else:
            deltat = 0.
        ue = Expression((initial_fields[ui]),
                        # element=VV[ui].ufl_element(),
                        degree=4,
                        t=t - deltat, nu=nu)
        NS_expressions["bc_%s" % ui] = ue

    bcs = dict((ui, []) for ui in sys_comp)
    bc0 = DirichletBC(V, NS_expressions["bc_u0"], boundary, 1)
    bc1 = DirichletBC(V, NS_expressions["bc_u1"], boundary, 1)

    bcs['u0'] = [bc0]
    bcs['u1'] = [bc1]
    bcs['p'] = []

    return bcs


def initialize(q_, q_1, q_2, VV, t, nu, dt, initial_fields, **NS_namespace):
    """Initialize solution.

    Use t=dt/2 for pressure since pressure is computed in between timesteps.

    """
    for ui in q_:
        if 'IPCS' in NS_parameters['solver']:
            deltat = dt / 2. if ui == 'p' else 0.
        else:
            deltat = 0.
        vv = interpolate(Expression((initial_fields[ui]),
                                    degree=4,
                                    t=t - deltat, nu=nu), VV[ui])
        q_[ui].vector()[:] = vv.vector()[:]
        if not ui == 'p':
            q_1[ui].vector()[:] = vv.vector()[:]
            deltat = -dt
            vv = interpolate(Expression((initial_fields[ui]),
                                        degree=4,
                                        t=t - deltat, nu=nu), VV[ui])
            q_2[ui].vector()[:] = vv.vector()[:]
    q_1['p'].vector()[:] = q_['p'].vector()[:]


class Walls(UserExpression):
    def __init__(self, t, coor, A, L, T_G, x_hat_map, counter_max, **kwargs):
        self.t = t
        self.map = x_hat_map
        self.max = counter_max
        self.coor = coor
        self.A = A
        self.L = L
        self.T_G = T_G
        self.counter = -1
        super().__init__(**kwargs)

    def update(self, t):
        self.t = t

    def eval(self, values, _, **kwargs):
        self.counter += 1
        if self.coor == 0:
            values[:] = 2 * np.pi / self.T_G * self.A * np.cos(2 * np.pi * self.t / self.T_G) * sin(
                2 * np.pi * (self.map[self.counter][1] + self.L / 2) / self.L)
        else:
            values[:] = 2 * np.pi / self.T_G * self.A * np.cos(2 * np.pi * self.t / self.T_G) * sin(
                2 * np.pi * (self.map[self.counter][0] + self.L / 2) / self.L)
        if self.counter == self.max:
            self.counter = -1


def pre_solve_hook(V, mesh, t, nu, L, w_, T_G, A0, newfolder, velocity_degree, u_components, boundary, exact_fields,
                   **NS_namespace):
    # Create exact solution
    ue_x = Expression(exact_fields['u0'], nu=nu, t=t,
                      degree=4)
    ue_y = Expression(exact_fields['u1'], nu=nu, t=t, degree=4)
    pe = Expression(exact_fields['p'], nu=nu, t=t, degree=4)

    # Visualization files
    viz_p, viz_u, viz_ue = get_visualization_writers(newfolder, ["Pressure", "Velocity", "Exact"])

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")
    ue_vec = Function(VV, name="Velocity")

    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    counter_max, x_hat_map = get_coordinate_map(V, boundary, w_, 1)
    walls_x = Walls(0, 0, A0, L, T_G, x_hat_map, counter_max, element=V.ufl_element())
    walls_y = Walls(0, 1, A0, L, T_G, x_hat_map, counter_max, element=V.ufl_element())

    NS_expressions["bc_w0"] = walls_x
    NS_expressions["bc_w1"] = walls_y

    # Mesh velocity conditions
    bc_mesh = dict((ui, []) for ui in u_components)
    bc0 = DirichletBC(V, walls_x, boundary, 1)
    bc1 = DirichletBC(V, walls_y, boundary, 1)

    bc_mesh["u0"] = [bc0]
    bc_mesh["u1"] = [bc1]
    return dict(viz_p=viz_p, viz_ue=viz_ue, ue_vec=ue_vec, viz_u=viz_u, u_vec=u_vec, dof_map=dof_map,
                bc_mesh=bc_mesh, coordinates=coordinates, ue_x=ue_x, ue_y=ue_y, pe=pe)


def update_boundary_conditions(t, dt, NS_expressions, **NS_namespace):
    for key, value in NS_expressions.items():
        if "bc" in key:
            if 'IPCS' in NS_parameters['solver'] and 'p' in key:
                deltat_ = dt / 2.
            else:
                deltat_ = 0.
            if "w" in key:
                value.update(t - deltat_)
            else:
                value.t = t - deltat_


def temporal_hook(q_, t, nu, VV, dt, u_vec, ue_vec, p_, viz_u, viz_p, viz_ue, initial_fields, tstep,
                  sys_comp, compute_error, total_error, ue_x, ue_y, pe, **NS_namespace):
    """Function called at end of timestep.

    Plot solution and compute error by comparing to analytical solution.
    Remember pressure is computed in between timesteps.

    """
    # Save solution
    assign(u_vec.sub(0), q_["u0"])
    assign(u_vec.sub(1), q_["u1"])

    viz_u.write(u_vec, t)
    viz_p.write(p_, t)

    ue_x.t = t
    ue_y.t = t
    pe.t = t - dt / 2
    uxx = interpolate(ue_x, VV['u0'])
    uyy = interpolate(ue_y, VV['u1'])
    pp = interpolate(pe, VV['p'])
    ues = []
    if tstep % compute_error == 0:
        err = {}
        for i, ui in enumerate(sys_comp):
            if 'IPCS' in NS_parameters['solver']:
                deltat_ = dt / 2. if ui == 'p' else 0.
            else:
                deltat_ = 0.
            ue = Expression((initial_fields[ui]),
                            element=VV[ui].ufl_element(),
                            t=t - deltat_, nu=nu)

            if ui == 'u0':
                ue = uxx
                assign(ue_vec.sub(0), ue)
            elif ui == 'u1':
                ue = uyy
                assign(ue_vec.sub(1), ue)
            elif ui == 'p':
                ue = pp
                pp.rename("p", 'p')
                viz_ue.write(pp, t)

            if "u" in ui:
                ues.append(ue)
            uen = norm(ue.vector())
            ue.vector().axpy(-1, q_[ui].vector())
            error = norm(ue.vector()) / uen
            err[ui] = "{0:2.6e}".format(norm(ue.vector()) / uen)
            total_error[i] += error * dt


def theend_hook(newfolder, mesh, q_, t, dt, nu, VV, sys_comp, total_error, initial_fields, **NS_namespace):
    final_error = np.zeros(len(sys_comp))
    for i, ui in enumerate(sys_comp):
        if 'IPCS' in NS_parameters['solver'] and ui == "p":
            deltat = dt / 2.
        else:
            deltat = 0.
        ue = Expression((initial_fields[ui]),
                        degree=4,
                        # element=VV[ui].ufl_element(),
                        t=t - deltat, nu=nu)
        ue = interpolate(ue, VV[ui])
        final_error[i] = errornorm(q_[ui], ue)

    hmin = mesh.hmin()
    if MPI.rank(MPI.comm_world) == 0:
        print("hmin = {}".format(hmin))
    s0 = "Total Error:"
    s1 = "Final Error:"
    for i, ui in enumerate(sys_comp):
        s0 += " {0:}={1:2.6e}".format(ui, total_error[i])
        s1 += " {0:}={1:2.6e}".format(ui, final_error[i])

    if MPI.rank(MPI.comm_world) == 0:
        print(s0)
        print(s1)

    err_u = final_error[0]
    err_ux = final_error[0]
    err_uy = final_error[1]
    err_p = final_error[2]
    err_p_h = final_error[2]

    # Write errors to file
    error_array = np.asarray([round(t, 4), err_u, err_ux, err_uy, err_p, err_p_h])
    error_path = path.join(newfolder, "error_mms.txt")
    with open(error_path, 'a') as filename:
        filename.write(
            f"{error_array[0]} {error_array[1]} {error_array[2]} {error_array[3]} {error_array[4]} {error_array[5]}\n")
