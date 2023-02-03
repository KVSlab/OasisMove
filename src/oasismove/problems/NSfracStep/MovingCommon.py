from oasismove.problems.NSfracStep import *


def get_coordinate_map(V, boundary, w_, facet_id):
    box_counter = Counter(element=V.ufl_element())
    bc_tmp = DirichletBC(V, box_counter, boundary, facet_id)
    bc_tmp.apply(w_["u0"].vector())
    w_["u0"].vector().zero()
    counter_max = box_counter.counter
    x_hat_map = box_counter.get_map()

    return counter_max, x_hat_map


class Counter(UserExpression):
    def __init__(self, **kwargs):
        self.map = {}
        self.counter = -1
        super().__init__(**kwargs)

    def get_map(self):
        return self.map

    def eval(self, _, x):
        self.counter += 1
        if x.size == 2:
            self.map[self.counter] = [x[0], x[1]]
        if x.size == 3:
            self.map[self.counter] = [x[0], x[1], x[2]]


def get_visualization_files(newfolder):
    viz_u = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "velocity.xdmf"))
    viz_p = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "pressure.xdmf"))
    for viz in [viz_u, viz_p]:
        viz.parameters["rewrite_function_mesh"] = True
        viz.parameters["flush_output"] = True
        viz.parameters["functions_share_mesh"] = True
    return viz_p, viz_u


def get_visualization_writers(newfolder, list_of_quantities):
    # List of quantities = ['velocity',..]
    writers = []
    for item in list_of_quantities:
        viz = XDMFFile(MPI.comm_world, path.join(newfolder, "Solutions", "{}.xdmf".format(item)))
        viz.parameters["rewrite_function_mesh"] = True
        viz.parameters["flush_output"] = True
        viz.parameters["functions_share_mesh"] = True
        writers.append(viz)

    return writers
