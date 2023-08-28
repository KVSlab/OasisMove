from dolfin import MPI, XDMFFile, UnitSquareMesh, FunctionSpace, Function

comm = MPI.comm_world
u_path = "velocity.xdmf"

# Create function
mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
t = 0.0

viz = XDMFFile(comm, u_path)
# Code hangs here:
viz.write(u, t)
