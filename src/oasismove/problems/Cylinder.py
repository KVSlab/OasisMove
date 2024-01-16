__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2014-04-10"
__copyright__ = "Copyright (C) 2014 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from dolfin import AutoSubDomain, near



Um = 5
nu = 0.004/1000
R = 1
H = 2

# Specify boundary conditions
Inlet = AutoSubDomain(lambda x, on_bnd: on_bnd and near(x[2], 0))
Wall = AutoSubDomain(lambda x, on_bnd: on_bnd and near(x[1]*x[1] + x[0]*x[0], R**2, 1e-2))
Outlet = AutoSubDomain(lambda x, on_bnd: on_bnd and near(x[2], H))


# Overload post_import_problem to choose between the two cases
def post_import_problem(NS_parameters, mesh, commandline_kwargs, **NS_namespace):
    # If the mesh is a callable function, then create the mesh here.
    if callable(mesh):
        mesh = mesh(**NS_parameters)

    """ Choose case - case could be defined through command line."""
    NS_parameters.update(commandline_kwargs)

    NS_parameters.update(nu=nu, Um=Um, R=R, mesh=mesh)

    return NS_parameters
