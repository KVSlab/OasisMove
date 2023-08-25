from dolfin import MPI, XDMFFile, UnitSquareMesh, FunctionSpace, Function

comm = MPI.comm_world
u_path = "velocity.h5"
file_mode = "w"

# Create function
mesh = UnitSquareMesh(1, 1)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
t = 0.0

if MPI.rank(comm) == 0:
    print("Writing")

# viz_u = HDF5File(MPI.comm_world, u_path, file_mode=file_mode) # Hangs

viz = XDMFFile(MPI.comm_world, "u.xdmf")
viz.write(u, t)  # Hangs

if MPI.rank(comm) == 0:
    print("Done")
