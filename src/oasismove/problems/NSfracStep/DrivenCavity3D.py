# Written by Mikael Mortensen <mikaem@math.uio.no> (2013)
# Edited by Henrik Kjeldsberg <henrik.kjeldsberg@live.no> (2023)
import pickle
from os import getcwd

from oasismove.problems.NSfracStep import *


# Override some problem specific parameters
def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(getcwd(), restart_folder)

        f = open(path.join(path.dirname(path.abspath(__file__)), restart_folder,
                           'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)

    else:
        NS_parameters.update(
            # Mesh parameters
            Nx=15,
            Ny=15,
            Nz=15,
            # Fluid parameters
            nu=0.01,
            # Simulation parameters
            T=1.0,
            dt=0.01,
            folder="results_driven_cavity_3d",
            # Oasis parameters
            testing=False,
            max_iter=1,
            dynamic_mesh=False,
            save_solution_frequency=5e10,
            save_step=2,
            checkpoint=50,
            print_intermediate_info=100,
            velocity_degree=1,
            pressure_degree=1,
            use_krylov_solvers=True)

    NS_expressions.update(dict(constrained_domain=PeriodicDomain()))


# Create a mesh
def mesh(Nx, Ny, Nz, **params):
    mesh = BoxMesh(Point(0, 0, 0), Point(1, 1, 1), Nx, Ny, Nz)

    return mesh


class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two slave edges
        return bool(near(x[2], 0) and on_boundary)

    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2] - 1.0


def create_bcs(V, **NS_namespace):
    # Specify boundary conditions
    noslip = "std::abs(x[0]*x[1]*(1-x[0]))<1e-8"
    top = "std::abs(x[1]-1) < 1e-8"

    bc0 = DirichletBC(V, 0, noslip)
    bc00 = DirichletBC(V, 1, top)
    bc01 = DirichletBC(V, 0, top)

    return dict(u0=[bc00, bc0],
                u1=[bc01, bc0],
                u2=[bc01, bc0],
                p=[])


def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_2:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
        [bc.apply(x_2[ui]) for bc in bcs[ui]]
