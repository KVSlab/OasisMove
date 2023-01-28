__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2014-04-10"
__copyright__ = "Copyright (C) 2014 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from dolfin import UnitSquareMesh


# Create a mesh
def mesh(Nx=50, Ny=50, **params):
    m = UnitSquareMesh(Nx, Ny)
    return m


noslip = "std::abs(x[0]*x[1]*(1-x[0]))<1e-8"
top = "std::abs(x[1]-1) < 1e-8"
bottom = "std::abs(x[1]) < 1e-8"
