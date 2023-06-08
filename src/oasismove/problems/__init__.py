# Written by Mikael Mortensen <mikaem@math.uio.no> (2013)
# Edited by Henrik Kjeldsberg <henrik.kjeldsberg@live.no> (2023)

import subprocess
from collections import defaultdict
from os import getpid, path
from numpy import array, maximum, zeros, savetxt
import numpy as np
from dolfin import *


def getMemoryUsage(rss=True):
    mypid = str(getpid())
    rss = "rss" if rss else "vsz"
    process = subprocess.Popen(['ps', '-o', rss, mypid],
                               stdout=subprocess.PIPE)
    out, _ = process.communicate()
    mymemory = out.split()[1]
    return eval(mymemory) / 1024


parameters["linear_algebra_backend"] = "PETSc"
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3"

# Default parameters for all solvers
NS_parameters = dict(
    nu=0.01,  # Kinematic viscosity
    folder='results',  # Relative folder for storing results
    velocity_degree=2,  # default velocity degree
    pressure_degree=1  # default pressure degree
)

NS_expressions = {}

constrained_domain = None

# To solve for scalars provide a list like ['scalar1', 'scalar2']
scalar_components = []

# With diffusivities given as a Schmidt number defined by:
#   Schmidt = nu / D (= momentum diffusivity / mass diffusivity)
Schmidt = defaultdict(lambda: 1.)
Schmidt_T = defaultdict(lambda: 0.7)  # Turbulent Schmidt number (LES)

Scalar = defaultdict(lambda: dict(Schmidt=1.0,
                                  family="CG",
                                  degree=1))

# The following helper functions are available in dolfin
# They are redefined here for printing only on process 0.
RED = "\033[1;37;31m%s\033[0m"
BLUE = "\033[1;37;34m%s\033[0m"
GREEN = "\033[1;37;32m%s\033[0m"


def info_blue(s, check=True):
    if MPI.rank(MPI.comm_world) == 0 and check:
        print(BLUE % s)


def info_green(s, check=True):
    if MPI.rank(MPI.comm_world) == 0 and check:
        print(GREEN % s)


def info_red(s, check=True):
    if MPI.rank(MPI.comm_world) == 0 and check:
        print(RED % s)


class OasisTimer(Timer):
    def __init__(self, task, verbose=False):
        Timer.__init__(self, task)
        info_blue(task, verbose)


class OasisMemoryUsage:
    def __init__(self, s):
        self.memory = 0
        self.memory_vm = 0
        self(s)

    def __call__(self, s, verbose=True):
        self.prev = self.memory
        self.prev_vm = self.memory_vm
        self.memory = MPI.sum(MPI.comm_world, getMemoryUsage())
        self.memory_vm = MPI.sum(MPI.comm_world, getMemoryUsage(False))
        if MPI.rank(MPI.comm_world) == 0 and verbose:
            info_blue('{0:26s}  {1:10d} MB {2:10d} MB {3:10d} MB {4:10d} MB'
                      .format(s,
                              int(self.memory - self.prev),
                              int(self.memory),
                              int(self.memory_vm - self.prev_vm),
                              int(self.memory_vm))
                      )


# Print memory use up til now
initial_memory_use = getMemoryUsage()
oasis_memory = OasisMemoryUsage('Start')


# Convenience functions
def strain(u):
    return 0.5 * (grad(u) + grad(u).T)


def omega(u):
    return 0.5 * (grad(u) - grad(u).T)


def Omega(u):
    return inner(omega(u), omega(u))


def Strain(u):
    return inner(strain(u), strain(u))


def QC(u):
    return Omega(u) - Strain(u)


def compute_flow_quantities(u, L, nu, mesh, t, tstep, dt, h, outlet_area=1, boundary=None, outlet_ids=[], inlet_ids=[],
                            id_wall=0, period=1.0, newfolder=None, dynamic_mesh=False, write_to_file=False):
    """
    Compute max and mean Reynolds number using velocities U_max & U_mean,
    kinematic viscosity (nu) and characteristic length (L)
    Also computes the CFL number, and fluxes through boundaries
    """

    # Compute and printCFL number
    DG = FunctionSpace(mesh, "DG", 0)
    U = project(sqrt(inner(u, u)), DG)

    U_mean = U.vector().get_local().mean()
    U_max = U.vector().get_local().max()

    re_mean = U_mean * L / nu
    re_max = U_max * L / nu
    flux_in = []
    flux_out = []
    flux_wall = []
    re_outlet = 0

    if boundary is not None:
        # Compute Reynolds number at outlet
        ds = Measure("ds", subdomain_data=boundary)
        n = FacetNormal(mesh)
        u_dot_n = dot(u, n)
        flux_wall.append(assemble(u_dot_n * ds(id_wall)))
        for inlet_id in inlet_ids:
            flux_in.append(assemble(u_dot_n * ds(inlet_id)))

        for outlet_id in outlet_ids:
            flux_out.append(assemble(u_dot_n * ds(outlet_id)))

        u_mean = flux_out[0] / outlet_area
        re_outlet = u_mean * L / nu

    # Compute Womersley number
    omega = 2 * np.pi / period
    D_outlet = np.sqrt(4 * outlet_area / np.pi)

    womersley_number_outlet = D_outlet * np.sqrt(omega / nu)

    # Recompute edge length if mesh has moved
    if dynamic_mesh:
        h = mesh.hmin()

    cfl = U.vector().get_local() * dt / h
    max_cfl = cfl.max()
    mean_cfl = cfl.mean()

    if MPI.rank(MPI.comm_world) == 0:
        info_green(
            'Time = {0:2.4e}, timestep = {1:6d}, max Reynolds number={2:2.3f}, mean Reynolds number={3:2.3f}, outlet Reynolds number={4:2.3f}, Womersley number={5:2.3f}, max velocity={6:2.3f}, mean velocity={7:2.3f}, max CFL={8:2.3f}, mean CFL={9:2.3f}'
                .format(t, tstep, re_max, re_mean, re_outlet, womersley_number_outlet, U_max, U_mean, max_cfl,
                        mean_cfl))

    if write_to_file:
        data = [t, tstep, re_max, re_mean, re_outlet, womersley_number_outlet, U_max, U_mean, max_cfl,
                mean_cfl] + flux_in + flux_out + flux_wall
        write_data_to_file(newfolder, data, "flow_metrics.txt")


def write_data_to_file(save_path, data, filename):
    data_path = path.join(save_path, filename)
    with open(data_path, "ab") as f:
        savetxt(f, data, fmt=" %.16f ", newline=' ')
        f.write(b'\n')


def print_mesh_information(mesh, dt=None, u_mean=None, dim=3):
    comm = MPI.comm_world
    local_xmin = mesh.coordinates()[:, 0].min()
    local_xmax = mesh.coordinates()[:, 0].max()
    local_ymin = mesh.coordinates()[:, 1].min()
    local_ymax = mesh.coordinates()[:, 1].max()
    if dim == 3:
        local_zmin = mesh.coordinates()[:, 2].min()
        local_zmax = mesh.coordinates()[:, 2].max()
    xmin = comm.gather(local_xmin, 0)
    xmax = comm.gather(local_xmax, 0)
    ymin = comm.gather(local_ymin, 0)
    ymax = comm.gather(local_ymax, 0)
    if dim == 3:
        zmin = comm.gather(local_zmin, 0)
        zmax = comm.gather(local_zmax, 0)

    local_hmin = mesh.hmin()
    local_num_cells = mesh.num_cells()
    local_num_edges = mesh.num_edges()
    local_num_faces = mesh.num_faces()
    local_num_facets = mesh.num_facets()
    local_num_vertices = mesh.num_vertices()
    h_min = comm.gather(local_hmin, 0)
    num_cells = comm.gather(local_num_cells, 0)
    num_edges = comm.gather(local_num_edges, 0)
    num_faces = comm.gather(local_num_faces, 0)
    num_facets = comm.gather(local_num_facets, 0)
    num_vertices = comm.gather(local_num_vertices, 0)
    volume = assemble(Constant(1) * dx(mesh))

    if MPI.rank(MPI.comm_world) == 0:
        print("=== Mesh information ===")
        print("X range: {} to {} (delta: {:.4f})".format(min(xmin), max(xmax), max(xmax) - min(xmin)))
        print("Y range: {} to {} (delta: {:.4f})".format(min(ymin), max(ymax), max(ymax) - min(ymin)))
        if dim == 3:
            print("Z range: {} to {} (delta: {:.4f})".format(min(zmin), max(zmax), max(zmax) - min(zmin)))
        print("Number of cells: {}".format(sum(num_cells)))
        print("Number of cells per processor: {}".format(int(np.mean(num_cells))))
        print("Number of edges: {}".format(sum(num_edges)))
        print("Number of faces: {}".format(sum(num_faces)))
        print("Number of facets: {}".format(sum(num_facets)))
        print("Number of vertices: {}".format(sum(num_vertices)))
        print("Volume: {:.4f}".format(volume))
        print("Number of cells per volume: {:.4f}".format(sum(num_cells) / volume))
        if u_mean is not None and dt is not None:
            print("CFL number: {:.4f}".format(u_mean * dt / min(h_min)))


def recursive_update(dst, src):
    """Update dict dst with items from src deeply ("deep update")."""
    for key, val in src.items():
        if key in dst and isinstance(val, dict) and isinstance(dst[key], dict):
            dst[key] = recursive_update(dst[key], val)
        else:
            dst[key] = val
    return dst


class OasisXDMFFile(XDMFFile, object):
    def __init__(self, comm, filename):
        XDMFFile.__init__(self, comm, filename)


def add_function_to_tstepfiles(function, newfolder, tstepfiles, tstep):
    name = function.name()
    tstepfolder = path.join(newfolder, "Timeseries")
    tstepfiles[name] = OasisXDMFFile(MPI.comm_world,
                                     path.join(tstepfolder,
                                               '{}_from_tstep_{}.xdmf'.format(name, tstep)))
    tstepfiles[name].function = function
    tstepfiles[name].parameters["rewrite_function_mesh"] = False


def body_force(mesh, **NS_namespace):
    """Specify body force"""
    return Constant((0,) * mesh.geometry().dim())


def initialize(**NS_namespace):
    """Initialize solution."""
    pass


def create_bcs(sys_comp, **NS_namespace):
    """Return dictionary of Dirichlet boundary conditions."""
    return dict((ui, []) for ui in sys_comp)


def scalar_hook(**NS_namespace):
    """Called prior to scalar solve."""
    pass


def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    return dict((ci, Constant(0)) for ci in scalar_components)


def pre_solve_hook(**NS_namespace):
    """Called just prior to entering time-loop. Must return a dictionary."""
    return {}


def theend_hook(**NS_namespace):
    """Called at the very end."""
    pass


def problem_parameters(**NS_namespace):
    """Updates problem specific parameters, and handles restart"""
    pass


def post_import_problem(NS_parameters, mesh, commandline_kwargs, NS_expressions, **NS_namespace):
    """Called after importing from problem."""
    dynamic_mesh = NS_parameters['dynamic_mesh']
    restart_folder = NS_parameters['restart_folder']
    # Update NS_parameters with all parameters modified through command line
    for key, val in commandline_kwargs.items():
        if isinstance(val, dict):
            NS_parameters[key].update(val)
        else:
            NS_parameters[key] = val

    # If the mesh is a callable function, then create the mesh here.
    if callable(mesh):
        mesh = mesh(**NS_parameters)

    # Load mesh if restarting simulation
    if restart_folder is not None:
        # Get mesh information
        mesh = Mesh()
        filename = path.join(restart_folder, 'mesh.h5')
        with HDF5File(MPI.comm_world, filename, "r") as mesh_file:
            mesh_file.read(mesh, "mesh", False)

    assert (isinstance(mesh, Mesh))

    # Returned dictionary to be updated in the NS namespace
    d = dict(mesh=mesh)
    d.update(NS_parameters)
    d.update(NS_expressions)
    return d


def u_dot_n(u, n):
    return (dot(u, n) - abs(dot(u, n))) / 2


def compute_volume(mesh, t, newfolder):
    volume = assemble(Constant(1.0) * dx(mesh))
    data = [t, volume]
    write_data_to_file(newfolder, data, "volume.txt")
