import pickle
from os import getcwd

from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_coordinate_map, get_visualization_writers


def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    """
    Problem file for running CFD simulation for the MovingVortex problem inspired by the problem by Fehn et al.[1],
    resembling the 2D Taylor-Green vortex. The problem solves the N-S equations in the absence of body forces, with a
    manufactured velocity solution. The mesh velocity is also described by an analytic displacement field,
    describing the oscillatory boundary movement. The movement is mainly controlled by the amplitude A0,
    and period length of the mesh motion T_G.

    [1] Fehn, N., Heinz, J., Wall, W. A., & Kronbichler, M. (2021). High-order arbitrary Lagrangian–Eulerian
    discontinuous Galerkin methods for the incompressible Navier–Stokes equations.
    Journal of Computational Physics, 430, 110040.
    """
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(getcwd(), restart_folder)
        f = open(path.join(path.dirname(path.abspath(__file__)), restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)
    else:
        # Read final time from command line
        T = float(commandline_kwargs.get("T", 1))
        NS_parameters.update(
            # Problem specific parameters
            A0=0.08,  # Amplitude
            T_G=4 * T,  # Period time
            L=1.0,  # Dimension of domain
            Nx=20,  # Resolution in x-direction
            Ny=20,  # Resolution in y-direction
            # Fluid parameters
            nu=0.025,  # Kinematic viscosity
            # Simulation parameters
            T=T,  # Simulation time
            dt=0.05,  # Time step
            folder="results_moving_vortex",
            # Oasis parameters
            max_iter=2,
            dynamic_mesh=True,
            save_step=1,
            save_solution_frequency=5,
            checkpoint=500,
            print_intermediate_info=100,
            compute_error=1,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True
        )

        NS_parameters['krylov_solvers'] = {'monitor_convergence': False,
                                           'report': False,
                                           'relative_tolerance': 1e-8,
                                           'absolute_tolerance': 1e-8}

    # Define analytical and initial state expressions
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


def mesh(Nx, Ny, L, dt, **params):
    # Define mesh
    mesh = RectangleMesh(Point(-L / 2, -L / 2), Point(L / 2, L / 2), Nx, Ny)

    print_mesh_information(mesh, dt, u_mean=1.0, dim=2)
    return mesh


def pre_boundary_condition(mesh, **NS_namespace):
    # Mark geometry
    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    inlet = AutoSubDomain(lambda x, b: b)

    inlet.mark(boundary, 1)
    return dict(boundary=boundary)


def create_bcs(V, t, dt, nu, sys_comp, boundary, initial_fields, **NS_namespace):
    for i, ui in enumerate(sys_comp):
        if 'IPCS' in NS_parameters['solver']:
            deltat = dt / 2. if ui == 'p' else 0.
        else:
            deltat = 0.
        ue = Expression((initial_fields[ui]),
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
        """
        (User)Expression class for describing the wall motion.
        """
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
        N = (self.coor + 1) % 2
        values[:] = 2 * np.pi / self.T_G * self.A * np.cos(2 * np.pi * self.t / self.T_G) * sin(
            2 * np.pi * (self.map[self.counter][N] + self.L / 2) / self.L)
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
    viz_p, viz_u, viz_ue = get_visualization_writers(newfolder, ["pressure", "velocity", "velocity_exact"])

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


def temporal_hook(u_, q_, t, nu, VV, dt, u_vec, ue_vec, p_, viz_u, viz_p, viz_ue, initial_fields, tstep,
                  save_solution_frequency, sys_comp, compute_error, total_error, ue_x, ue_y, pe, testing,
                  **NS_namespace):
    """Function called at end of timestep.

    Plot solution and compute error by comparing to analytical solution.
    Remember pressure is computed in between timesteps.

    """
    # Save solution
    if not testing and tstep % save_solution_frequency == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])

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

            if "u" in ui:
                ues.append(ue)
            uen = norm(ue.vector())
            ue.vector().axpy(-1, q_[ui].vector())
            error = norm(ue.vector()) / uen
            err[ui] = "{0:2.6e}".format(norm(ue.vector()) / uen)
            total_error[i] += error * dt

        viz_ue.write(ue_vec, t)


def theend_hook(q_, t, dt, nu, VV, sys_comp, total_error, initial_fields, **NS_namespace):
    final_error = np.zeros(len(sys_comp))
    for i, ui in enumerate(sys_comp):
        if 'IPCS' in NS_parameters['solver'] and ui == "p":
            deltat = dt / 2.
        else:
            deltat = 0.
        ue = Expression((initial_fields[ui]),
                        element=VV[ui].ufl_element(),
                        t=t - deltat, nu=nu)
        ue = interpolate(ue, VV[ui])
        final_error[i] = errornorm(q_[ui], ue)

    s0 = "Total Error:"
    s1 = "Final Error:"
    for i, ui in enumerate(sys_comp):
        s0 += " {0:}={1:2.6e}".format(ui, total_error[i])
        s1 += " {0:}={1:2.6e}".format(ui, final_error[i])

    if MPI.rank(MPI.comm_world) == 0:
        print(s0)
        print(s1)
