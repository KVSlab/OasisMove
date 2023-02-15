# flake8: noqa
import matplotlib.pyplot as plt
from IPython import embed
from oasis.problems.NSfracStep.MovingCommon import get_visualization_writers

from .LagrangianParticleTracking import LagrangianParticles
from .particle_generators import RandomRectangle
from ..NSfracStep import *


# Override some problem specific parameters
def problem_parameters(NS_parameters, scalar_components, Schmidt, NS_expressions, **NS_namespace):
    Re = 100.
    U_in = 1
    L = 0.5
    nu = U_in * L / Re
    T = 30
    NS_parameters.update(dict(
        U_in=U_in,
        nu=nu,
        T=T,
        dt=0.1,
        Re=Re,
        Nx=40,
        Ny=40,
        folder="results_moving_backstep",
        max_iter=2,
        plot_interval=1,
        print_intermediate_info=1e10,
        velocity_degree=1,
        save_step=1,
        cfl_step=10,  # Determine how often CFL step is computed, printed, and saved
        compute_re_step=10,  # Determine how often Reynolds number (Re) is computed, printed and saved
        compute_flux_step=10,  # Determine how often Reynolds number (Re) is computed, printed and saved
        use_krylov_solvers=True))

    scalar_components += ["blood"]
    Schmidt["blood"] = 100


def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    print(scalar_components)
    return dict((ci, Constant(1)) for ci in scalar_components)


# Create a mesh here
class Submesh(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] <= 2 - DOLFIN_EPS and x[1] <= 0.5 - DOLFIN_EPS


def mesh(Ny, **params):
    Nx = 5 * Ny
    mesh_ = RectangleMesh(Point(0, 0), Point(5, 1), Nx, Ny)
    # return mesh_
    subm = Submesh()
    mf1 = MeshFunction("size_t", mesh_, 2)
    mf1.set_all(0)
    subm.mark(mf1, 1)
    return SubMesh(mesh_, mf1, 0)


def pre_boundary_condition(mesh, **NS_namespace):
    walls = AutoSubDomain(lambda x, b: b and x[0] <= 5.0 + DOLFIN_EPS)
    inlet = AutoSubDomain(lambda x, b: b and near(x[0], 0.0))
    outlet = AutoSubDomain(lambda x, b: b and x[0] >= 5.0 - DOLFIN_EPS)
    walls_top = AutoSubDomain(lambda x, b: b and x[1] > 0.75)

    boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary.set_all(0)

    walls.mark(boundary, 1)
    walls_top.mark(boundary, 2)
    inlet.mark(boundary, 3)
    outlet.mark(boundary, 4)

    return dict(boundary=boundary)


class Wall(UserExpression):
    def __init__(self, **kwargs):
        self.t = 0
        self.A = 1 / 2
        self.f = np.pi / 10
        self.fx = np.pi / 3
        self.U0 = 0
        super().__init__(**kwargs)

    def update(self, t):
        self.t = t
        self.F = np.sin(self.f * self.t)
        self.dF = self.f * self.A * np.cos(self.f * self.t)

        self.U0 = self.F / np.abs(self.F) * self.dF

    def eval(self, value, x):
        sigma = np.sqrt(2) / 3
        x0 = 5 / 2
        scale = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(- (x[0] - x0) ** 2 / (2 * sigma ** 2))
        value[:] = scale * self.U0


def create_bcs(V, Q, sys_comp, U_in, boundary, dynamic_mesh, **NS_namespace):
    walls_x = Constant(0)
    if dynamic_mesh:
        walls_y = Wall(element=V.ufl_element())
    else:
        walls_y = Constant(0)

    NS_expressions["walls_x"] = walls_x
    NS_expressions["walls_y"] = walls_y

    bcs = dict((ui, []) for ui in sys_comp)
    # Bottom wall
    bc0 = DirichletBC(V, 0., boundary, 1)
    # Top wall
    bct_x = DirichletBC(V, walls_x, boundary, 2)
    bct_y = DirichletBC(V, walls_y, boundary, 2)
    # Inlet
    bcu_x = DirichletBC(V, U_in, boundary, 3)
    bcu_y = DirichletBC(V, 0, boundary, 3)
    # Outlet
    bcp = DirichletBC(Q, Constant(0.0), boundary, 4)
    bcb = DirichletBC(V, Constant(0.0), boundary, 3)
    bcs['u0'] = [bc0, bcu_x, bct_x]
    bcs['u1'] = [bc0, bcu_y, bct_y]
    bcs['p'] = [bcp]
    bcs['blood'] = [bcb]

    return bcs


def pre_solve_hook(mesh, newfolder, T, dt, save_step, nu, u_components, boundary, x_, V, velocity_degree,
                   **NS_namespace):
    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    up_vec = Function(Vv, name="up")

    quantities = ["blood", "velocity"]
    viz_b, viz_u = get_visualization_writers(newfolder, quantities)

    viz_lpt = []
    tstep = int(T / dt)
    n_tsteps = int(tstep / save_step)
    for i in range(n_tsteps):
        save_path = "particles_%03d.xdmf" % int(dt * (i + 1) * 10)
        viz = XDMFFile(MPI.comm_world, path.join(newfolder, "Particles", save_path))
        viz.parameters["rewrite_function_mesh"] = True
        viz.parameters["flush_output"] = True
        viz.parameters["functions_share_mesh"] = True
        viz_lpt.append(viz)

    # Compute edge length
    h = mesh.hmin()

    # Extract dof map and coordinates
    VV = VectorFunctionSpace(mesh, "CG", velocity_degree)
    u_vec = Function(VV, name="Velocity")

    # Extract dof map and coordinates
    coordinates = mesh.coordinates()
    if velocity_degree == 1:
        dof_map = vertex_to_dof_map(V)
    else:
        dof_map = V.dofmap().entity_dofs(mesh, 0)

    # Mesh velocity boundary conditions
    noslip = Constant(0.0)

    bcu_walls = DirichletBC(V, noslip, boundary, 1)

    bcu_walls_x = DirichletBC(V, noslip, boundary, 2)
    bcu_walls_y = DirichletBC(V, NS_expressions['walls_y'], boundary, 2)

    bc_in_x = DirichletBC(V, noslip, boundary, 3)
    bc_in_y = DirichletBC(V, noslip, boundary, 3)

    bc_out_x = DirichletBC(V, noslip, boundary, 4)
    bc_out_y = DirichletBC(V, noslip, boundary, 4)

    bc_mesh = dict((ui, []) for ui in u_components)

    rigid_bc_x = [bcu_walls_x, bcu_walls]
    rigid_bc_y = [bcu_walls_y, bcu_walls]
    bc_mesh["u0"] = [bc_in_x, bc_out_x] + rigid_bc_x
    bc_mesh["u1"] = [bc_in_y, bc_out_y] + rigid_bc_y

    # Setup particle tracking
    y0 = 0.5
    particle_positions = RandomRectangle(Point(0.01, y0), Point(0.1, 1)).generate([25, 25], method='full')

    # Create particles path in case running in parallel
    save_particles_path = path.join(newfolder, "particles.npy")
    lp = LagrangianParticles(Vv, save_particles_path)

    lp.add_particles(particle_positions, properties_d={"time": np.zeros(len(particle_positions))})

    Vv = VectorFunctionSpace(mesh, "CG", velocity_degree)
    w_vec = Function(Vv, name="u")

    return dict(viz_u=viz_u, viz_lpt=viz_lpt, w_vec=w_vec, u_vec=u_vec,
                bc_mesh=bc_mesh, dof_map=dof_map, coordinates=coordinates, lp=lp, h=h,
                viz_b=viz_b, up_vec=up_vec)


def update_boundary_conditions(t, dynamic_mesh, NS_expressions, **NS_namespace):
    # Update time
    for key, value in NS_expressions.items():
        if "walls_y" in key and dynamic_mesh:
            value.update(t)


def temporal_hook(tstep, newfolder, q_, viz_b, mesh, viz_lpt, t, viz_u, save_step, u_, u_vec, lp, dt, **NS_namespace):
    assign(u_vec.sub(0), u_[0])
    assign(u_vec.sub(1), u_[1])

    lp.step(u_vec, step=tstep, dt=dt)

    if tstep % save_step == 0:
        File(path.join(newfolder, "mesh_%04d.pvd" % int(10 * t))) << mesh
        viz_u.write(u_vec, int(10 * t))
        viz_b.write(q_['blood'], int(10 * t))
        lp.save_particles(viz_lpt[tstep - 1])


def theend_hook(lp, **NS_namespace):
    size = MPI.comm_world.Get_size()
    if size > 1:
        print("-- Saving particles")
        with open(lp.save_particles_path, 'wb') as f:
            np.save(f, np.array(lp.particles_history))
