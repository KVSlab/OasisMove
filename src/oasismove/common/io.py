# Written by Mikael Mortensen <mikaem@math.uio.no> (2013)
# Edited by Henrik Kjeldsberg <henrik.kjeldsberg@live.no> (2023)


import glob
import pickle
import time
from os import makedirs, listdir, remove, system, path
from xml.etree import ElementTree as ET

from dolfin import MPI, XDMFFile, HDF5File, FunctionSpace, Function, interpolate, MeshFunction

from oasismove.problems import info_red

__all__ = ["create_initial_folders", "save_solution", "save_tstep_solution_xdmf",
           "save_checkpoint_solution_xdmf", "check_if_kill", "check_if_reset_statistics",
           "init_from_restart", "merge_visualization_files", "merge_xml_files"]


def create_initial_folders(folder, restart_folder, sys_comp, tstep, info_red,
                           scalar_components, output_timeseries_as_vector,
                           **NS_namespace):
    """Create necessary folders."""
    info_red("Creating initial folders")
    # To avoid writing over old data create a new folder for each run
    if MPI.rank(MPI.comm_world) == 0:
        try:
            makedirs(folder)
        except OSError:
            pass

    MPI.barrier(MPI.comm_world)
    newfolder = path.join(folder, 'data')
    if restart_folder:
        newfolder = path.join(newfolder, restart_folder.split('/')[-2])
    else:
        if not path.exists(newfolder):
            newfolder = path.join(newfolder, '1')
        else:
            previous = [f for f in listdir(newfolder) if not f.startswith('.')]
            previous = max(map(eval, previous)) if previous else 0
            newfolder = path.join(newfolder, str(previous + 1))

    MPI.barrier(MPI.comm_world)
    if MPI.rank(MPI.comm_world) == 0:
        if not restart_folder:
            makedirs(path.join(newfolder, "Timeseries"))
            makedirs(path.join(newfolder, "Checkpoint"))

    tstepfolder = path.join(newfolder, "Timeseries")
    tstepfiles = {}
    comps = sys_comp
    if output_timeseries_as_vector:
        comps = ['p', 'u'] + scalar_components

    for ui in comps:
        tstepfiles[ui] = XDMFFile(MPI.comm_world, path.join(
            tstepfolder, ui + '_from_tstep_{}.xdmf'.format(tstep)))
        tstepfiles[ui].parameters["functions_share_mesh"] = True
        tstepfiles[ui].parameters["rewrite_function_mesh"] = True
        tstepfiles[ui].parameters["flush_output"] = True

    return newfolder, tstepfiles


def save_solution(tstep, t, q_, q_1, w_, d_, folder, newfolder, save_step, checkpoint, NS_parameters, tstepfiles,
                  u_, u_components, output_timeseries_as_vector, mesh, AssignedVectorFunction, killtime, total_timer,
                  **NS_namespace):
    """Called at end of timestep. Check for kill and save solution if required."""
    NS_parameters.update(t=t, tstep=tstep)
    if tstep % save_step == 0:
        save_tstep_solution_xdmf(tstep, q_, u_, newfolder, tstepfiles, output_timeseries_as_vector,
                                 AssignedVectorFunction, NS_parameters)

    pauseoasis = check_if_pause(folder)
    while pauseoasis:
        time.sleep(5)
        pauseoasis = check_if_pause(folder)

    killoasis = check_if_kill(folder, killtime, total_timer)
    if tstep % checkpoint == 0 or killoasis:
        save_checkpoint_solution_xdmf(q_, q_1, w_, d_, newfolder, u_components, mesh, NS_parameters)

    return killoasis


def save_tstep_solution_xdmf(tstep, q_, u_, newfolder, tstepfiles, output_timeseries_as_vector, AssignedVectorFunction,
                             NS_parameters):
    """Store solution on current timestep to XDMF file."""
    timefolder = path.join(newfolder, 'Timeseries')
    if output_timeseries_as_vector:
        # project or store velocity to vector function space
        for comp, tstepfile in tstepfiles.items():
            if comp == "u":
                # Create vector function and assigners
                uv = AssignedVectorFunction(u_)

                # Assign solution to vector
                uv()

                # Store solution vector
                tstepfile.write(uv, float(tstep))

            elif comp in q_:
                tstepfile.write(q_[comp], float(tstep))

            else:
                tstepfile.write(tstepfile.function, float(tstep))

    else:
        for comp, tstepfile in tstepfiles.items():
            tstepfile << (q_[comp], float(tstep))

    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(path.join(timefolder, "params.dat")):
            f = open(path.join(timefolder, 'params.dat'), 'wb')
            pickle.dump(NS_parameters, f)


def save_checkpoint_solution_xdmf(q_, q_1, w_, d_, newfolder, u_components, mesh, NS_parameters):
    """
    Overwrite solution in Checkpoint folder
    """
    checkpointfolder = path.join(newfolder, "Checkpoint")
    NS_parameters["num_processes"] = MPI.size(MPI.comm_world)
    if MPI.rank(MPI.comm_world) == 0:
        if path.exists(path.join(checkpointfolder, "params.dat")):
            system('cp {0} {1}'.format(path.join(checkpointfolder, "params.dat"),
                                       path.join(checkpointfolder, "params_old.dat")))
        f = open(path.join(checkpointfolder, "params.dat"), 'wb')
        pickle.dump(NS_parameters, f)

    if MPI.rank(MPI.comm_world) == 0 and path.exists(path.join(checkpointfolder, "params_old.dat")):
        system('rm {0}'.format(path.join(checkpointfolder, "params_old.dat")))

    # Store velocity and pressure solution
    MPI.barrier(MPI.comm_world)
    for ui in q_:
        checkpoint_path = path.join(checkpointfolder, ui + '.xdmf')
        with XDMFFile(MPI.comm_world, checkpoint_path) as f:
            f.write_checkpoint(q_[ui], '/current')
            if ui in u_components:
                f.write_checkpoint(q_1[ui], '/previous', append=True)
        MPI.barrier(MPI.comm_world)

    MPI.barrier(MPI.comm_world)
    # Store mesh velocity and deformation solution
    if w_ is not None:
        for ui in w_:
            checkpoint_path = path.join(checkpointfolder, ui.replace("u", "w") + '.xdmf')
            with XDMFFile(MPI.comm_world, checkpoint_path) as f:
                f.write_checkpoint(w_[ui], '/current')

            checkpoint_path = path.join(checkpointfolder, ui.replace("u", "d") + '.xdmf')
            with XDMFFile(MPI.comm_world, checkpoint_path) as f:
                f.write_checkpoint(d_[ui], '/current')
            MPI.barrier(MPI.comm_world)

    # Store mesh and boundary
    MPI.barrier(MPI.comm_world)
    mesh_path = path.join(checkpointfolder, 'mesh.h5')
    with HDF5File(MPI.comm_world, mesh_path, 'w') as f:
        f.write(mesh, 'mesh')
        boundary = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())
        f.write(boundary, 'boundary')


def check_if_kill(folder, killtime, total_timer):
    """Check if user has put a file named killoasis in folder or if given killtime has been reached."""
    found = 0
    if 'killoasis' in listdir(folder):
        found = 1
    collective = MPI.sum(MPI.comm_world, found)
    if collective > 0:
        if MPI.rank(MPI.comm_world) == 0:
            remove(path.join(folder, 'killoasis'))
            info_red('killoasis Found! Stopping simulations cleanly...')
        return True
    else:
        elapsed_time = float(total_timer.elapsed()[0])
        if killtime is not None and killtime <= elapsed_time:
            if MPI.rank(MPI.comm_world) == 0:
                info_red('Given killtime reached! Stopping simulations cleanly...')
            return True
        else:
            return False


def check_if_pause(folder):
    """Check if user has put a file named pauseoasis in folder."""
    found = 0
    if 'pauseoasis' in listdir(folder):
        found = 1
    collective = MPI.sum(MPI.comm_world, found)
    if collective > 0:
        if MPI.rank(MPI.comm_world) == 0:
            info_red('pauseoasis Found! Simulations paused. Remove '
                     + path.join(folder, 'pauseoasis') + ' to resume simulations...')
        return True
    else:
        return False


def check_if_reset_statistics(folder):
    """Check if user has put a file named resetoasis in folder."""
    found = 0
    if 'resetoasis' in listdir(folder):
        found = 1
    collective = MPI.sum(MPI.comm_world, found)
    if collective > 0:
        if MPI.rank(MPI.comm_world) == 0:
            remove(path.join(folder, 'resetoasis'))
            info_red('resetoasis Found!')
        return True
    else:
        return False


def init_from_restart(restart_folder, sys_comp, uc_comp, u_components, q_, q_1, q_2, w_, d_, tstep, velocity_degree,
                      previous_velocity_degree, mesh, constrained_domain, V, Q, **NS_namespace):
    """Initialize solution from checkpoint files """
    if restart_folder:
        if MPI.rank(MPI.comm_world) == 0:
            info_red('Restarting from checkpoint at time step {}'.format(tstep))
        q_prev = q_
        q_2_prev = q_2
        if previous_velocity_degree != velocity_degree:
            # Create dictionaries for the solutions at previous and different element degree
            V_prev = FunctionSpace(mesh, 'CG', previous_velocity_degree, constrained_domain=constrained_domain)
            VV_prev = dict((ui, V_prev) for ui in uc_comp)
            VV_prev['p'] = Q

            q_prev = dict((ui, Function(VV_prev[ui], name=ui)) for ui in sys_comp)
            q_2_prev = dict((ui, Function(V_prev, name=ui + "_2")) for ui in u_components)

        # Load velocity and pressure
        for ui in sys_comp:
            checkpoint_path = path.join(restart_folder, ui + '.xdmf')
            with XDMFFile(MPI.comm_world, checkpoint_path) as f:
                # Interpolate
                read_and_interpolate_solution(f, V, previous_velocity_degree, q_, q_prev, ui, velocity_degree,
                                              "/current")
                if ui in uc_comp:
                    q_1[ui].vector().zero()
                    q_1[ui].vector().axpy(1., q_[ui].vector())
                    q_1[ui].vector().apply('insert')
                    if ui in u_components:
                        # Interpolate
                        read_and_interpolate_solution(f, V, previous_velocity_degree, q_2, q_2_prev, ui,
                                                      velocity_degree, "/previous")
        # Load mesh velocity and deformation
        if w_ is not None:
            for ui in u_components:
                checkpoint_path = path.join(restart_folder, ui.replace("u", "w") + '.xdmf')
                with XDMFFile(MPI.comm_world, checkpoint_path) as f:
                    f.read_checkpoint(w_[ui], "/current")

                checkpoint_path = path.join(restart_folder, ui.replace("u", "d") + '.xdmf')
                with XDMFFile(MPI.comm_world, checkpoint_path) as f:
                    f.read_checkpoint(d_[ui], "/current")


def read_and_interpolate_solution(f, V, previous_velocity_degree, q_, q_prev, ui, velocity_degree, name):
    """
    Interpolate solution to higher element order or read directly into existing function space
    """
    if previous_velocity_degree != velocity_degree and ui != 'p':
        f.read_checkpoint(q_prev[ui], name)
        q_prev_proj = interpolate(q_prev[ui], V)
        q_[ui].vector().zero()
        q_[ui].vector().axpy(1., q_prev_proj.vector())
        q_[ui].vector().apply('insert')
    else:
        f.read_checkpoint(q_[ui], name)


def merge_visualization_files(newfolder, **namesapce):
    timefolder = path.join(newfolder, 'Timeseries')
    # Gather files
    xdmf_files = list(glob.glob(path.join(timefolder, "*.xdmf")))
    xdmf_velocity = [f for f in xdmf_files if "u_from_tstep" in f.__str__()]
    xdmf_pressure = [f for f in xdmf_files if "p_from_tstep" in f.__str__()]

    # Merge files
    for files in [xdmf_velocity, xdmf_pressure]:
        if len(files) > 1:
            merge_xml_files(files)


def merge_xml_files(files):
    # Get first timestep and trees
    first_timesteps = []
    trees = []
    MPI.barrier(MPI.comm_world)
    for f in files:
        trees.append(ET.parse(f))
        root = trees[-1].getroot()
        first_timesteps.append(float(root[0][0][0][2].attrib["Value"]))
        MPI.barrier(MPI.comm_world)

    # Index valued sort (bypass numpy dependency)
    first_timestep_sorted = sorted(first_timesteps)
    indexes = [first_timesteps.index(i) for i in first_timestep_sorted]

    # Get last timestep of first tree
    base_tree = trees[indexes[0]]
    last_node = base_tree.getroot()[0][0][-1]
    ind = 1 if len(list(last_node)) == 3 else 2
    last_timestep = float(last_node[ind].attrib["Value"])

    # Append
    for index in indexes[1:]:
        tree = trees[index]
        for node in list(tree.getroot()[0][0]):
            ind = 1 if len(list(node)) == 3 else 2
            if last_timestep < float(node[ind].attrib["Value"]):
                base_tree.getroot()[0][0].append(node)
                last_timestep = float(node[ind].attrib["Value"])

    # Seperate xdmf files
    new_file = [f for f in files if "_0" in f]
    old_files = [f for f in files if "_" in f and f not in new_file]

    # Write new xdmf file
    base_tree.write(new_file[0], xml_declaration=True)

    # Delete xdmf file
    if MPI.rank(MPI.comm_world) == 0:
        [remove(f) for f in old_files]
