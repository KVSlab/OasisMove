from dolfin import MPI, XDMFFile, UnitSquareMesh, FunctionSpace, Function, HDF5File

comm = MPI.comm_world
u_path = "velocity.xdmf"

# Create function
mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
t = 0.0

# Code hangs here:

viz = XDMFFile(comm, u_path)
viz.write(u, t)  # Hangs

if MPI.rank(comm) == 0:
    print("Done")
