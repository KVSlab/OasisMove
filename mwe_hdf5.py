from dolfin import MPI, HDF5File

comm = MPI.comm_world
u_path = "velocity.h5"
file_mode = "w"

# Code hangs here:
viz_u = HDF5File(comm, u_path, file_mode=file_mode)
