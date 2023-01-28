__author__ = 'Anna Haley <ahaley@mie.utoronto.ca>'
__date__ = '2020-05-04'
__copyright__ = 'Copyright (C) 2020 ' + __author__
__license__ = 'GNU Lesser GPL version 3 or any later version'

from dolfin import (SpatialCoordinate, exp, sqrt)

from oasismove.problems.NSfracStep.MovingAtriumCommon import get_distance_to_plane

__all__ = ['nn_setup', 'nn_update']

import numpy as np
from ufl import conditional, gt


def nn_setup(mesh, p_MV, p_FE, FE_rad, nu, CG1Function, nu_nn_krylov_solver, **NS_namespace):
    """
    Set up for solving Modified-Cross non-Newtonian model.
    """
    # Set up Modified Cross form
    n = p_FE - p_MV
    n /= np.linalg.norm(n)
    x = SpatialCoordinate(mesh)
    p = np.array([x[0], x[1], x[2]])

    dist_flowext = get_distance_to_plane(n, p_FE, p)

    c1 = 0.2
    c2 = dist_flowext
    nu_dist = 20 * nu * exp(-c1 * c2) / (1 + exp(-c1 * c2))

    # Limit the viscosity change to a sphere around the MV
    rad = sqrt((p[0] - p_FE[0]) ** 2 + (p[1] - p_FE[1]) ** 2 + (p[2] - p_FE[2]) ** 2)
    dist_laa = conditional(gt(rad, FE_rad * 1.2), 0, 1)
    nu_nn_form = nu + nu_dist * dist_laa
    nunn_ = CG1Function(nu_nn_form, mesh, method=nu_nn_krylov_solver, bounded=True, name="nu_nn")

    # Call nunn_ only once
    nunn_()
    return dict(nunn_=nunn_, bcs_nu_nn=[])


def nn_update(nunn_, **NS_namespace):
    """Compute nunn_"""
    pass  # nunn_()
