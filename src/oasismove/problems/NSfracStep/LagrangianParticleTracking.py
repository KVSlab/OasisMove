# flake8: noqa
# __authors__ = ('Mikael Mortensen <mikaem@math.uio.no>',
#                'Miroslav Kuchta <mirok@math.uio.no>',
#                'Henrik Kjeldsberg <henriakj@simula.no>')
# __date__ = '2022-20-02'
# __copyright__ = 'Copyright (C) 2011' + __authors__
# __license__  = 'GNU Lesser GPL version 3 or any later version'
"""
This module contains functionality for Lagrangian tracking of particles with
DOLFIN
"""

import dolfin as df
import numpy as np
import copy
from mpi4py import MPI as pyMPI
from collections import defaultdict

# Disable printing
__DEBUG__ = False

comm = pyMPI.COMM_WORLD

# collision tests return this value or -1 if there is no collision
__UINT32_MAX__ = np.iinfo('uint32').max


class Particle:
    __slots__ = ['position', 'properties']
    'Lagrangian particle with position and some other passive properties.'

    def __init__(self, x):
        self.position = x
        self.properties = {}

    def send(self, dest):
        """Send particle to dest."""
        comm.Send(self.position, dest=dest)
        comm.send(self.properties, dest=dest)

    def recv(self, source):
        """Receive info of a new particle sent from source."""
        comm.Recv(self.position, source=source)
        self.properties = comm.recv(source=source)


class CellWithParticles(df.Cell):
    """Dolfin cell with list of particles that it contains."""

    def __init__(self, mesh, cell_id, particle):
        # Initialize parent -- create Cell with id on mesh
        df.Cell.__init__(self, mesh, cell_id)
        # Make an empty list of particles that I carry
        self.particles = []
        self += particle

    def __add__(self, particle):
        """Add single particle to cell."""
        assert isinstance(particle, (Particle, np.ndarray))
        if isinstance(particle, Particle):
            self.particles.append(particle)
            return self
        else:
            return self.__add__(Particle(particle))

    def __len__(self):
        """Number of particles in cell."""
        return len(self.particles)


class CellParticleMap(dict):
    """Dictionary of cells with particles."""

    def __add__(self, ins):
        """
        Add ins to map:
            ins is either (mesh, cell_id, particle) or
                          (mesh, cell_id, particle, particle_properties)
        """
        assert isinstance(ins, tuple) and len(ins) in (3, 4)
        # If the cell_id is in map add the particle
        if ins[1] in self:
            self[ins[1]] += ins[2]
        # Other wise create new cell
        else:
            self[ins[1]] = CellWithParticles(ins[0], ins[1], ins[2])
        # With particle_properties, update properties of the last added particle
        if len(ins) == 4:
            self[ins[1]].particles[-1].properties.update(ins[3])

        return self

    def pop(self, cell_id, i):
        """Remove i-th particle from the list of particles in cell with cell_id."""
        # Note that we don't check for cell_id being a key or cell containg
        # at least i particles.
        particle = self[cell_id].particles.pop(i)

        # If the cell is empty remove it from map
        if len(self[cell_id]) == 0:
            del self[cell_id]

        return particle

    def total_number_of_particles(self):
        """Total number of particles in all cells of the map."""
        return sum(map(len, iter(self.values())))


class LagrangianParticles:
    """Particles moved by the velocity field in V."""

    def __init__(self, V, particles_path):
        self.__debug = __DEBUG__

        self.save_particles_path = particles_path
        self.V = V
        self.mesh = V.mesh()
        self.mesh.init(2, 2)  # Cell-cell connectivity for neighbors of cell
        self.tree = self.mesh.bounding_box_tree()  # Tree for section comput.
        self.particles_history = []

        # Allocate some variables used to look up the velocity is computed as U_i*basis_i where i is the dimension of
        # element function space, U are coefficients and basis_i are element
        # function space basis functions. For interpolation in cell it is
        # advantageous to compute the restriction once for cell and only
        # update basis_i(x) depending on x, i.e. particle where we make
        # interpolation. This update mounts to computing the basis matrix
        self.dim = self.mesh.topology().dim()

        self.element = V.dolfin_element()
        self.num_tensor_entries = 1

        rank = len(V.ufl_element().value_shape())

        for i in range(rank):
            self.num_tensor_entries *= self.element.value_dimension(i)
        # For VectorFunctionSpace CG1 this is 3
        self.coefficients = np.zeros(self.element.space_dimension())
        # For VectorFunctionSpace CG1 this is 3x3
        self.basis_matrix = np.zeros((self.element.space_dimension(),
                                      self.num_tensor_entries))

        # Allocate a dictionary to hold all particles
        self.particle_map = CellParticleMap()

        # Allocate some MPI stuff
        self.num_processes = comm.Get_size()
        self.my_rank = comm.Get_rank()
        self.all_processes = list(range(self.num_processes))
        self.other_processes = list(range(self.num_processes))
        self.other_processes.remove(self.my_rank)
        self.my_escaped_particles = np.zeros(1, dtype='I')
        self.tot_escaped_particles = np.zeros(self.num_processes, dtype='I')
        # Dummy particle for receiving/sending at [0, 0, ...]
        self.particle0 = Particle(np.zeros(self.mesh.geometry().dim()))

    def __iter__(self):
        """Iterate over all particles."""
        for cwp in self.particle_map.values():
            for particle in cwp.particles:
                yield particle

    def add_particles(self, list_of_particles, properties_d=None):
        """Add particles and search for their home on all processors.
           Note that list_of_particles must be same on all processes. Further
           every len(properties[property]) must equal len(list_of_particles).
        """
        my_found, all_found = self.insert_particles(list_of_particles, properties_d)
        # All particles must be found on some process
        comm.Reduce(my_found, all_found, root=0)

        if self.my_rank == 0:
            missing = np.where(all_found == 0)[0]
            n_missing = len(missing)
            new_props = {}
            new_points = []
            for ID in missing:
                new_props = {"time": np.zeros(len(missing))}
                pt = list_of_particles[ID]
                #rng = np.random.uniform(-0.249, 0.249, size=1)[0]
                #pt.position = np.array([0.01, 0.75 + rng])
                rng = np.random.uniform(-0.49, 0.49, size=1)[0]
                pt.position = np.array([0.01, 0.5 + rng])
                new_points.append(pt)

            self.insert_particles(new_points, properties_d=new_props)

            # Print particle info
            if self.__debug:
                for i in missing:
                    print('Missing', list_of_particles[i].position)

                n_duplicit = len(np.where(all_found > 1)[0])
                print('There are %d duplicit particles' % n_duplicit)

    def insert_particles(self, list_of_particles, properties_d):
        """Add particles and search for their home on all processors.
           Note that list_of_particles must be same on all processes. Further
           every len(properties[property]) must equal len(list_of_particles).
        """
        if properties_d is not None:
            n = len(list_of_particles)
            assert all(len(sub_list) == n
                       for sub_list in properties_d.values())
            # Dictionary that will be used to feed properties of single
            # particles
            properties = list(properties_d.keys())
            particle_properties = dict((key, 0) for key in properties)

            has_properties = True
        else:
            has_properties = False

        pmap = self.particle_map
        my_found = np.zeros(len(list_of_particles), 'I')
        all_found = np.zeros(len(list_of_particles), 'I')

        for i, particle in enumerate(list_of_particles[::-1]):
            c = self.locate(particle)
            if not (c == -1 or c == __UINT32_MAX__):
                my_found[i] = True
                if not has_properties:
                    pmap += self.mesh, c, particle
                else:
                    # Get values of properties for this particle
                    for key in properties:
                        particle_properties[key] = properties_d[key][i]
                    pmap += self.mesh, c, particle, particle_properties

        return my_found, all_found

    def step(self, u, step, dt):
        """Move particles by forward Euler x += u*dt"""

        # Update bounding box in case mesh moved
        self.tree = self.mesh.bounding_box_tree()
        self.tree.build(self.mesh)

        start = df.Timer('shift')
        for cwp in self.particle_map.values():
            # Restrict once per cell
            self.coefficients[:] = u.restrict(self.element, cwp)
            for particle in cwp.particles:
                x = particle.position  # Compute velocity at position x
                cell_vtx = cwp.get_vertex_coordinates()
                cell_orient = cwp.orientation()

                for i, row in enumerate(self.basis_matrix):
                    row[:] = self.element.evaluate_basis(i, x, cell_vtx, cell_orient)
                x[:] = x[:] + dt * np.dot(self.coefficients, self.basis_matrix)[:]
                particle.properties['time'] = particle.properties['time'] + dt

        # Recompute the map
        stop_shift = start.stop()
        start = df.Timer('relocate')
        self.relocate()
        stop_reloc = start.stop()
        # We return computation time per process
        return stop_shift, stop_reloc

    def relocate(self):
        # Relocate particles on cells and processors
        p_map = self.particle_map
        # Map such that map[old_cell] = [(new_cell, particle_id), ...]
        # Ie new destination of particles formerly in old_cell
        new_cell_map = defaultdict(list)
        for cwp in p_map.values():
            for i, particle in enumerate(cwp.particles):
                point = df.Point(*particle.position)
                # Search only if particle moved outside original cell
                if not cwp.contains(point):
                    found = False
                    # Check neighbor cells
                    for neighbor in df.cells(cwp):
                        if neighbor.contains(point):
                            new_cell_id = neighbor.index()
                            found = True
                            break
                    # Do a completely new search if not found by now
                    if not found:
                        new_cell_id = self.locate(particle)
                    # Record to map
                    new_cell_map[cwp.index()].append((new_cell_id, i))

        # Rebuild locally the particles that end up on the process. Some
        # have cell_id == -1, i.e. they are on other process
        list_of_escaped_particles = []
        for old_cell_id, new_data in new_cell_map.items():
            # We iterate in reverse because normal order would remove some
            # particle this shifts the whole list!
            for (new_cell_id, i) in sorted(new_data,
                                           key=lambda t: t[1],
                                           reverse=True):
                particle = p_map.pop(old_cell_id, i)
                if new_cell_id == -1 or new_cell_id == __UINT32_MAX__:
                    list_of_escaped_particles.append(particle)
                else:
                    p_map += self.mesh, new_cell_id, particle

        # Create a list of how many particles escapes from each processor
        self.my_escaped_particles[0] = len(list_of_escaped_particles)
        # Make all processes aware of the number of escapees
        comm.Allgather(self.my_escaped_particles, self.tot_escaped_particles)

        # Send particles to root
        if self.my_rank != 0:
            for particle in list_of_escaped_particles:
                particle.send(0)

        # Receive the particles escaping from other processors
        if self.my_rank == 0:
            for proc in self.other_processes:
                for i in range(self.tot_escaped_particles[proc]):
                    self.particle0.recv(proc)
                    list_of_escaped_particles.append(copy.deepcopy(self.particle0))

        # Put all travelling particles on all processes, then perform new search
        travelling_particles = comm.bcast(list_of_escaped_particles, root=0)
        self.add_particles(travelling_particles)

    def total_number_of_particles(self):
        """Return number of particles in total and on process."""
        num_p = self.particle_map.total_number_of_particles()
        tot_p = comm.allreduce(num_p)
        return tot_p, num_p

    def locate(self, particle):
        """Find mesh cell that contains particle."""
        assert isinstance(particle, (Particle, np.ndarray))
        if isinstance(particle, Particle):
            # Convert particle to point
            point = df.Point(*particle.position)
            return self.tree.compute_first_entity_collision(point)
        else:
            return self.locate(Particle(particle))

    def save_particles(self, viz_lpt):
        """Scatter plot of all particles on process 0"""
        p_map = self.particle_map
        all_particles = np.zeros(self.num_processes, dtype='I')
        my_particles = p_map.total_number_of_particles()
        # Root learns about count of particles on all processes
        comm.Gather(np.array([my_particles], 'I'), all_particles, root=0)

        # Slaves should send to master
        if self.my_rank > 0:
            for cwp in p_map.values():
                for p in cwp.particles:
                    p.send(0)
        else:
            # Receive on master
            received = defaultdict(list)
            received_property = defaultdict(list)
            received[0] = [copy.copy(p.position)
                           for cwp in p_map.values()
                           for p in cwp.particles]
            received_property[0] = [copy.copy(p.properties['time'])
                                    for cwp in p_map.values()
                                    for p in cwp.particles]
            for proc in self.other_processes:
                # Receive all_particles[proc]
                for j in range(all_particles[proc]):
                    self.particle0.recv(proc)
                    received[proc].append(copy.copy(self.particle0.position))
                    received_property[proc].append(copy.copy(self.particle0.properties['time']))

            self.save_position(viz_lpt, received, received_property)

    def save_position(self, viz_lpt, particles, properties):
        property_idx = "time"
        size = comm.Get_size()
        all_particles = []
        all_properties = []
        for received, props in zip(particles.values(), properties.values()):
            all_particles += list(received)
            all_properties += list(props)
        # Workaround to save in parallel. Saving list does not seem to work in parallel.
        if size > 1:
            # Dump to numpy array
            particles_data = [[p[0], p[1], pr] for p, pr in zip(all_particles, all_properties)]
            self.particles_history.append(particles_data)
        else:
            all_particles = []
            for received in particles.values():
                all_particles += list(received)

            points_list = [df.Point(p[0], p[1]) for p in all_particles]
            viz_lpt.write(points_list, all_properties)

    def bar(self, fig):
        """Bar plot of particle distribution."""
        ax = fig.gca()

        p_map = self.particle_map
        all_particles = np.zeros(self.num_processes, dtype='I')
        my_particles = p_map.total_number_of_particles()
        # Root learns about count of particles on all processes
        comm.Gather(np.array([my_particles], 'I'), all_particles, root=0)

        if self.my_rank == 0 and self.num_processes > 1:
            ax.bar(np.array(self.all_processes) - 0.25, all_particles, 0.5)
            ax.set_xlabel('proc')
            ax.set_ylabel('number of particles')
            ax.set_xlim(-0.25, max(self.all_processes) + 0.25)
            return np.sum(all_particles)
        else:
            return None
