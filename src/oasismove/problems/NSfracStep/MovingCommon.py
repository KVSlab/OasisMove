from os import path

from dolfin import DirichletBC
from dolfin import UserExpression, XDMFFile, MPI


def get_coordinate_map(V, boundary, w_, facet_id):
    """
     Create a mapping of coordinates for a given function space and boundary.

     Args:
         V (dolfin.function.Functionspace): A FunctionSpace object representing the field on which the map is created.
         boundary (dolfin.cpp.mesh.MeshFunction): A MeshFunction defining with marked boundaries for the mesh.
         w_ (dict): A dict containing the mesh velocity at the current time step.
         facet_id (int): An integer representing the id of the boundary facets.

     Returns:
         counter_max (int): An integer count of the number of elements in the field `V`.
         x_hat_map (dict): A dict mapping coordinates to elements in the field `V`.
     """
    box_counter = Counter(element=V.ufl_element())
    bc_tmp = DirichletBC(V, box_counter, boundary, facet_id)
    bc_tmp.apply(w_["u0"].vector())
    w_["u0"].vector().zero()

    return box_counter.counter, box_counter.get_map()


class Counter(UserExpression):
    """
    A UserExpression designed for evaluating and counting mesh nodes,
    and associating each node with a unique integer, thus keeping
    a reference to the original mesh.

    Attributes:
        map (dict): A dictionary that maps a unique integer to a node coordinate.
        counter (int): An integer that counts the number of evaluated nodes.
    """

    def __init__(self, **kwargs):
        self.map = {}
        self.counter = -1
        super().__init__(**kwargs)

    def get_map(self):
        """ Returns the map that associates each unique integer with a node coordinate. """
        return self.map

    def eval(self, _, x):
        """
        Evaluates the node coordinate and increments the counter. The node coordinate is stored in the map with the
        counter value as the key.

        Args:
            x (numpy.ndarray): An array representing the node coordinate.
        """
        self.counter += 1
        self.map[self.counter] = list(x)


def get_visualization_writers(newfolder, list_of_quantities):
    """
    Create a list of XDMFFile writers for different quantities.

    Args:
      newfolder (str): The path of the folder where the XDMF files will be saved.
      list_of_quantities (list): A list of string names representing different quantities to be written,
      such as 'velocity' or 'pressure'.

    Returns:
      writers (list): A list of XDMFFile objects for writing each quantity to a separate file.
    """
    writers = []
    for item in list_of_quantities:
        viz = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "{}.xdmf".format(item)))
        viz.parameters["rewrite_function_mesh"] = True
        viz.parameters["flush_output"] = True
        viz.parameters["functions_share_mesh"] = True
        writers.append(viz)

    return writers
