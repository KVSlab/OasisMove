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
            T=0.01,
            dt=0.01,
            checkpoint=1,
            save_step=1,
            plot_interval=10,
            velocity_degree=2,
            print_intermediate_info=100,
            use_krylov_solvers=True,
            mesh_path=commandline_kwargs["mesh_path"],
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


def create_bcs(V, Q, Um, R, **NS_namespace):
    inlet = Expression(
        f"Um*(R*R-(x[0]*x[0]+x[1]*x[1]))/(R*R)", R=R, Um=Um, degree=2)
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


def temporal_hook(mesh, V, newfolder, tstepfiles, tstep, ds, u_,
                   AssignedVectorFunction, **NS_namespace):
    global_uv = AssignedVectorFunction(u_, name='Velocity')
    global_uv()
    xdmffile.write_checkpoint(global_uv, "u", tstep, XDMFFile.Encoding.HDF5)
    #omega = Function(V, name='omega')
    # Store omega each save_step
    #add_function_to_tstepfiles(omega, newfolder, tstepfiles, tstep)
    
    
    # ff = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    # Inlet.mark(ff, 1)
    # Outlet.mark(ff, 2)
    # Wall.mark(ff, 3)
    #breakpoint()
    #n = FacetNormal(mesh)
    #ds = ds[ff]
    return {}
    #return dict(uv=uv, omega=omega, ds=ds, ff=ff, n=n)


#def temporal_hook(q_, u_, tstep, V, uv, p_, plot_interval, omega, ds,
 #                 save_step, mesh, nu,, n, **NS_namespace):
# if tstep % plot_interval == 0:
#     uv()
#     plt.figure(1)
#     plot(uv, title='Velocity')
#     plt.figure(2)
#     plot(p_, title='Pressure')
#     plt.figure(3)
#     plot(q_['alfa'], title='alfa')
#     plt.show()

# R = VectorFunctionSpace(mesh, 'R', 0)
# c = TestFunction(R)
# tau = -p_ * Identity(2) + nu * (grad(u_) + grad(u_).T)
# forces = assemble(dot(dot(tau, n), c) * ds(1)).get_local() * 2 / Umean ** 2 / D

# print("Cd = {}, CL = {}".format(*forces))

# if tstep % save_step == 0:
#     try:
#         from fenicstools import StreamFunction
#         omega.assign(StreamFunction(u_, []))
#     except Exception:
#         omega.assign(project(curl(u_), V, solver_type='cg',
#                              bcs=[DirichletBC(V, 0, DomainBoundary())]))


#def theend_hook(q_, u_, p_, uv, mesh, ds, V, nu, Umean, D, **NS_namespace):
    # uv()
    # plot(uv, title='Velocity')
    # plot(p_, title='Pressure')
    # plot(q_['alfa'], title='alfa')
    # R = VectorFunctionSpace(mesh, 'R', 0)
    # c = TestFunction(R)
    # tau = -p_ * Identity(2) + nu * (grad(u_) + grad(u_).T)
    # ff = MeshFunction("size_t", mesh, mesh.ufl_cell().geometric_dimension() - 1)
    # Cyl.mark(ff, 1)
    # n = FacetNormal(mesh)
    # ds = ds[ff]
    # forces = assemble(dot(dot(tau, n), c) * ds(1)).get_local() * 2 / Umean ** 2 / D

    # print("Cd = {}, CL = {}".format(*forces))

    # try:
    #     from fenicstools import Probes
    #     from numpy import linspace, repeat, where, resize
    #     xx = linspace(0, L, 10000)
    #     x = resize(repeat(xx, 2), (10000, 2))
    #     x[:, 1] = 0.2
    #     probes = Probes(x.flatten(), V)
    #     probes(u_[0])
    #     nmax = where(probes.array() < 0)[0][-1]
    #     print("L = ", x[nmax, 0] - 0.25)
    #     print("dP = ", p_(Point(0.15, 0.2)) - p_(Point(0.25, 0.2)))

    # except Exception:
    #     pass

def theend_hook(**NS_namespace):
    xdmffile.close()