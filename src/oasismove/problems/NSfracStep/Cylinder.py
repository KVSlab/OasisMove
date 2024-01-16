from __future__ import print_function

__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2014-03-21"
__copyright__ = "Copyright (C) 2014 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

import pickle
from os import getcwd

import matplotlib.pyplot as plt
from oasismove.problems.NSfracStep import *
from oasismove.problems.Cylinder import *
from pathlib import Path

t=0
inlet = Expression(
        f"(t<1 ? t : 1)*Um*(R*R-(x[0]*x[0]+x[1]*x[1]))/(R*R)",R=R, Um=Um, t=t, degree=2)

def problem_parameters(commandline_kwargs, NS_parameters, scalar_components,
                       Schmidt, **NS_namespace):
    # Example: python NSfracstep.py [...] restart_folder="results/data/8/Checkpoint"
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(getcwd(), restart_folder)
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)
    else:
        # Override some problem specific parameters
        NS_parameters.update(
            T=2,
            dt=0.1,
            checkpoint=1,
            save_step=1,
            plot_interval=10,
            velocity_degree=2,
            print_intermediate_info=100,
            use_krylov_solvers=True,
            mesh_path=commandline_kwargs["mesh_path"],
            inlet=inlet
        )
        NS_parameters['krylov_solvers'].update(dict(monitor_convergence=True))
        NS_parameters['velocity_krylov_solver'].update(dict(preconditioner_type='jacobi',
                                                            solver_type='bicgstab'))

    scalar_components.append("alfa")
    Schmidt["alfa"] = 0.1


def mesh(mesh_path, dt, **NS_namespace):
    # Import mesh
    print(mesh_path)
    if Path(mesh_path).suffix == ".xml":
        mesh = Mesh(mesh_path)
    elif Path(mesh_path).suffix == ".xdmf":
        mesh = Mesh()
        with XDMFFile(mesh_path) as infile:
            infile.read(mesh)

    print_mesh_information(mesh, dt, dim=2)
    return mesh

def create_bcs(V, Q,**NS_namespace):

    ux = Expression("0.00", degree=2)
    uy = Expression("0.00", degree=2)
    bc_i0 = DirichletBC(V, 0, Inlet)
    bc_i1 = DirichletBC(V, 0, Inlet)
    bc_i2 = DirichletBC(V, inlet, Inlet)

    bc2 = DirichletBC(V, 0, Wall)
    bcp = DirichletBC(Q, 0, Outlet)

    return dict(u0=[bc_i0, bc2],
                u1=[bc_i1, bc2],
                u2=[bc_i2, bc2],
                p=[bcp],
                alfa=[])


def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:
        [bc.apply(x_2[ui]) for bc in bcs[ui]]

xdmffile = XDMFFile(MPI.comm_world, "results/u.xdmf")
first = [True]

def start_timestep_hook(inlet, t,  **NS_namespace):
    inlet.t =t

def temporal_hook(t, u_,
                   AssignedVectorFunction, **NS_namespace):
    global_uv = AssignedVectorFunction(u_, name='Velocity')
    global_uv()
    xdmffile.write_checkpoint(global_uv, "u", t, XDMFFile.Encoding.HDF5, append=not first[0])
    if first[0]:
        first[0] = False
 
def theend_hook(**NS_namespace):
    xdmffile.close()