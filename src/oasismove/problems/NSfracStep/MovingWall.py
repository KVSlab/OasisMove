from oasismove.problems.NSfracStep import *
from oasismove.problems.NSfracStep.MovingCommon import get_coordinate_map, get_visualization_writers

comm = MPI.comm_world


def problem_parameters(NS_parameters, **NS_namespace):
    """
    Problem file for running CFD simulation for the MovingWall problem inspired by the Wall-induced channel flow
    presented by Chnafa [1]. The problem considers flow in a long, straight, time-dependent domain, where the flow is
    induced by a moving wall at y=h(t). The motion is mainly controlled by the pulsation (sigma) and the ampliture
    of oscillations (eps), with an initial height of h0. The flow has an analytic solution for small Reynolds numbers,
    and may be used as a validation case.

    [1] Chnafa, C. (2014). Using image-based large-eddy simulations to investigate the intracardiac flow and its
    turbulent nature (Doctoral dissertation, Université Montpellier II-Sciences et Techniques du Languedoc).
    """
    NS_parameters.update(
        # Problem specifiic parameters
        h0=0.001,  # Initial height of wall
        eps=0.05,  # Amplitude
        sigma=2 * np.pi,  # Pulsation of the movement
        u_max=0.01,  # Maximum velocity
        # Mesh parameters
        scale=25,  # Ratio between computational length of x-direction and y-direction
        Nx=250,  # Resolution on x-direction
        Ny=20,  # Resolution on x-direction
        # Fluid parameters
        nu=8 * 10 ** (-7),  # Kinematic viscosity
        # Simulation properties
        T=1.0,  # End time
        dt=0.001,  # Time step
        folder="results_moving_wall",
        # Oasis parameters
        max_iter=2,
        dynamic_mesh=True,
        save_solution_frequency=1,
        checkpoint=500,
        print_intermediate_info=100,
        velocity_degree=1,
        pressure_degree=1,
        use_krylov_solvers=True,
        max_error=1e-8
    )


def mesh(eps, h0, scale, dt, Nx, Ny, u_max, **NS_namespace):
    X = scale * h0
    Y = h0 * (1 + eps)

    mesh = RectangleMesh(Point(0, 0), Point(X, Y), Nx, Ny)
    print_mesh_information(mesh, dt, u_max, dim=2)
    return mesh


def pre_boundary_condition(mesh, h0, eps, scale, **NS_namespace):
    # Mark geometry
    movingwall = AutoSubDomain(lambda x: near(x[1], h0 * (1 + eps)))
    outlet = AutoSubDomain(lambda x: near(x[0], scale * h0))
    leftwall = AutoSubDomain(lambda x: near(x[0], 0))
    bottomwall = AutoSubDomain(lambda x: near(x[1], 0))

    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    movingwall.mark(boundary, 1)
    leftwall.mark(boundary, 2)
    bottomwall.mark(boundary, 3)
    outlet.mark(boundary, 4)

    return dict(boundary=boundary)


def create_bcs(V, Q, w_, sigma, h0, eps, sys_comp, boundary, NS_expressions, **NS_namespace):
    counter_left, left_coordinate_map = get_coordinate_map(V, boundary, w_, 2)
    counter_right, right_coordinate_map = get_coordinate_map(V, boundary, w_, 4)

    # Create inlet expressions
    movingwall = MovingWall(0, sigma, h0, eps, element=V.ufl_element())
    leftwall = MovingSideWall(0, sigma, h0, eps, left_coordinate_map, counter_left, element=V.ufl_element())
    rightwall = MovingSideWall(0, sigma, h0, eps, right_coordinate_map, counter_right, element=V.ufl_element())
    noslip = Constant(0.0)

    # Store expressions
    NS_expressions["movingwall"] = movingwall
    NS_expressions["leftwall"] = leftwall
    NS_expressions["rightwall"] = rightwall
    NS_expressions["bottom"] = noslip

    # Velocity
    bcu_movingwall_x = DirichletBC(V, noslip, boundary, 1)
    bcu_movingwall_y = DirichletBC(V, movingwall, boundary, 1)

    bcu_leftwall_x = DirichletBC(V, noslip, boundary, 2)
    bcu_leftwall_y = DirichletBC(V, leftwall, boundary, 2)

    bcu_bottomwall = DirichletBC(V, noslip, boundary, 3)

    # Pressure
    bcp_out = DirichletBC(Q, Constant(0), boundary, 4)

    bcs = dict((ui, []) for ui in sys_comp)
    bcs['u0'] = [bcu_movingwall_x, bcu_leftwall_x]
    bcs['u1'] = [bcu_movingwall_y, bcu_bottomwall, bcu_leftwall_y]
    bcs["p"] = [bcp_out]

    return bcs


def pre_solve_hook(V, mesh, newfolder, velocity_degree, u_components, boundary, **NS_namespace):
    # Visualization files
    viz_p, viz_u = get_visualization_writers(newfolder, ['pressure', 'velocity'])

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Mesh velocity conditions
    noslip = Constant(0.0)

    bcu_walls_x = DirichletBC(V, noslip, boundary, 1)
    bcu_walls_y = DirichletBC(V, NS_expressions["movingwall"], boundary, 1)

    bc_left_wall_x = DirichletBC(V, noslip, boundary, 2)
    bc_left_wall_y = DirichletBC(V, NS_expressions["leftwall"], boundary, 2)

    bc_bottom_wall = DirichletBC(V, noslip, boundary, 3)

    bc_out_x = DirichletBC(V, noslip, boundary, 4)
    bc_out_y = DirichletBC(V, NS_expressions["rightwall"], boundary, 4)

    bc_mesh = dict((ui, []) for ui in u_components)

    rigid_bc_x = [bc_bottom_wall, bc_out_x]
    rigid_bc_y = [bc_bottom_wall, bc_out_y]
    bc_mesh["u0"] = [bcu_walls_x, bc_left_wall_x] + rigid_bc_x
    bc_mesh["u1"] = [bcu_walls_y, bc_left_wall_y] + rigid_bc_y

    return dict(viz_p=viz_p, viz_u=viz_u, u_vec=u_vec, bc_mesh=bc_mesh, dof_map=dof_map, coordinates=coordinates)


class MovingWall(UserExpression):
    def __init__(self, t, sigma, h0, eps, **kwargs):
        self.t = t
        self.sigma = sigma
        self.h0 = h0
        self.eps = eps
        self.value = 0
        super().__init__(**kwargs)

    def update(self, t):
        """
        Sinusidal wall motion from [1]

        [1] Chnafa, C. (2014). Using image-based large-eddy simulations to investigate the intracardiac flow and its
        turbulent nature (Doctoral dissertation, Université Montpellier II-Sciences et Techniques du Languedoc).
        """
        self.t = t
        self.value = - self.sigma * self.h0 * self.eps * np.sin(self.sigma * t)

    def eval(self, value, _):
        value[:] = self.value


class MovingSideWall(UserExpression):
    def __init__(self, t, sigma, h0, eps, x_hat_map, counter_max, **kwargs):
        self.t = t
        self.map = x_hat_map
        self.max = counter_max
        self.counter = -1
        self.sigma = sigma
        self.h0 = h0
        self.eps = eps
        self.value = 0
        super().__init__(**kwargs)

    def update(self, t):
        self.t = t

    def eval(self, values, _, **kwargs):
        self.counter += 1
        scaling = self.map[self.counter][1] / (self.h0 * (1 + self.eps))
        values[:] = - scaling * self.sigma * self.h0 * self.eps * np.sin(self.sigma * self.t)
        if self.counter == self.max:
            self.counter = -1


def update_boundary_conditions(t, NS_expressions, **NS_namespace):
    # Update time
    for key, value in NS_expressions.items():
        if "wall" in key:
            value.update(t)


def temporal_hook(mesh, dt, h0, eps, nu, t, tstep, save_solution_frequency, p_, u_, viz_u, viz_p, u_vec,
                  **NS_namespace):
    # Write solution at time t
    if tstep % save_solution_frequency == 0:
        assign(u_vec.sub(0), u_[0])
        assign(u_vec.sub(1), u_[1])

        viz_u.write(u_vec, t)
        viz_p.write(p_, t)

    # Compute mean velocity and Reynolds number at inlet
    if tstep % 10 == 0:
        h = mesh.hmin()
        L = h0 * (1 + eps)
        compute_flow_quantities(u_, L, nu, mesh, t, tstep, dt, h, outlet_area=1, boundary=None, outlet_ids=[],
                                inlet_ids=[], id_wall=0, period=1.0, newfolder=None, dynamic_mesh=False,
                                write_to_file=False)
