__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2014-04-08"
__copyright__ = "Copyright (C) 2014 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from ..DrivenCavity import *
from ..NSCoupled import *


# Override some problem specific parameters
def problem_parameters(NS_parameters, **NS_namespace):
    NS_parameters.update(nu=0.002, max_iter=100)


# Specify boundary conditions
def create_bcs(VQ, **NS_namespace):
    bc0 = DirichletBC(VQ.sub(0), (0, 0), noslip)
    bc1 = DirichletBC(VQ.sub(0), (1, 0), top)
    return dict(up=[bc0, bc1])


def theend_hook(u_, p_, up_, mesh, testing, **NS_namespace):
    up_.set_allow_extrapolation(True)
    if MPI.rank(MPI.comm_world) == 0 and testing:
        u_corner = up_((1, 1))[1]
        print("Velocity in corner = {0:2.6e}".format(u_corner))

    if not testing:
        plot(u_, title="Velocity")
        plot(p_, title="Pressure")

        try:
            import matplotlib.pyplot as plt
            from fenicstools import StreamFunction

            psi = StreamFunction(u_, [], mesh, use_strong_bc=True)
            plot(psi, title="Streamfunction")
            plt.show()
        except ImportError:
            pass
