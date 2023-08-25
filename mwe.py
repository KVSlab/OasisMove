from dolfin import MPI, HDF5File

comm = MPI.comm_world
u_path = "velocity.h5"
file_mode = "w"
for i in range(10):
    if MPI.rank(comm) == 0:
        print("Writing")

    MPI.barrier(comm)
    viz_u = HDF5File(MPI.comm_world, u_path, file_mode=file_mode)
    MPI.barrier(comm)

    if MPI.rank(comm) == 0:
        print("Done")
