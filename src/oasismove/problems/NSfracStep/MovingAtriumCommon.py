from os import path, makedirs

import numpy as np
from dolfin import *
from scipy.interpolate import splev, splrep


class Wall_motion(UserExpression):
    def __init__(self, t, motion, max_counter, direction, cycle, **kwargs):
        self.t = t
        self.motion = motion
        self.max_counter = max_counter
        self.counter = -1
        self.direction = direction
        self.cycle = cycle
        super().__init__(**kwargs)

    def eval(self, values, _):
        self.counter += 1
        # No motion for inlet/outlet
        values[:] = splev(self.t % self.cycle, self.motion[self.counter][self.direction], der=1)
        if self.counter == self.max_counter:
            self.counter = -1


class NuSpace(UserExpression):
    def __init__(self, nu_inf, p_mv, p_flowext, n, **kwargs):
        self.nu_inf = nu_inf
        self.p_mv = p_mv
        self.p_flowext = p_flowext
        self.n = n
        self.tol = 1E-1  # in [mm]
        self.dist_max = get_distance_between_points(self.p_mv, self.p_flowext) + self.tol
        super().__init__(**kwargs)

    def eval(self, value, x):
        p = np.array([x[0], x[1], x[2]])
        dist_mv = d = get_distance_to_plane(self.n, self.p_mv, p)
        dist_flowext = get_distance_to_plane(self.n, self.p_flowext, p)

        if dist_mv <= self.dist_max and dist_flowext <= self.dist_max:
            value[0] = self.nu_inf * (1 + 19 * d / self.dist_max)
        else:
            value[0] = self.nu_inf


def get_distance_between_points(p0, p1):
    return np.linalg.norm(np.array(p1) - np.array(p0))


def get_distance_to_plane(n, P, Q):
    return np.abs(n.dot(P - Q))


def get_file_paths(folder):
    # Create folder where data and solutions (velocity, mesh, pressure) is stored
    common_path = path.join(folder, "Solutions")
    if MPI.rank(MPI.comm_world) == 0:
        if not path.exists(common_path):
            makedirs(common_path)

    file_p = path.join(common_path, "p.h5")
    file_u = path.join(common_path, "u.h5")
    file_u_mean = path.join(common_path, "u_mean.h5")
    file_mesh = path.join(common_path, "mesh.h5")
    files = {"u": file_u, "u_mean": file_u_mean, "p": file_p, "mesh": file_mesh}

    return files


class Surface_counter(UserExpression):
    def __init__(self, points, cycle, flow_rate_type, **kwargs):
        self.motion = {}
        self.counter = -1
        self.points = points
        self.flow_rate_type = flow_rate_type
        self.time = np.linspace(0, cycle, self.points.shape[-1])
        super().__init__(**kwargs)

    def get_motion(self):
        return self.motion

    def eval(self, _, x):
        self.counter += 1
        index = np.argmin(np.sqrt(np.sum((self.points[:, :, 0] - np.array(x)) ** 2, axis=1)))
        # FIXME: Somehow predict good s to select?
        if "AF" in self.flow_rate_type:
            s = 0.01
        elif "SR" in self.flow_rate_type:
            s = 0.5
        elif "LA5" in self.flow_rate_type:
            s = 0.0025
        else:
            print("Not valid flow rate")
            exit()

        x_ = splrep(self.time, self.points[index, 0, :], s=s, per=True)
        y_ = splrep(self.time, self.points[index, 1, :], s=s, per=True)
        z_ = splrep(self.time, self.points[index, 2, :], s=s, per=True)

        self.motion[self.counter] = [x_, y_, z_]


class InletParabolic(UserExpression):
    def __init__(self, tstep, dt, N, n, Q_profile, center, R2, area, area_total, cycle, power, **kwargs):
        self.tstep = tstep
        self.dt = dt
        self.N = N
        self.Q_profile = Q_profile
        self.Q = 0
        self.center = center
        self.R2 = R2
        self.area = area
        self.area_total = area_total
        self.normal_component = n
        self.cycle = cycle
        self.power = power
        self.t = self.dt * self.tstep
        super().__init__(**kwargs)

    def update(self, tstep):
        self.tstep = tstep
        tstep = self.tstep % self.N
        self.t = self.dt * tstep
        self.Q = splev(self.t % self.cycle, self.Q_profile)

    def eval(self, value, x):
        factor = self.area / self.area_total
        # To satisfy flux, scale velocity. For degree n, scale with: C = (n+2)/n
        if self.power == 0:
            # For coarse mesh: boundary will affect flow rate
            # At boundary U=0, hence scale U0 with
            # (inner points - boundary points/2) / total points ^(-1) = 1 / 0.88
            # coarse_mesh_scale = 1 / 0.88
            scale = 1.0  # coarse_mesh_scale
        else:
            scale = (self.power + 2) / self.power

        U0 = scale * self.Q / self.area * factor

        x0 = self.center[0]
        x1 = self.center[1]
        x2 = self.center[2]
        R2 = self.R2
        power = self.power
        if power == 0:
            parabolic = U0
        else:
            R_power = np.sqrt(R2) ** power
            parabolic = U0 * (1 - ((x0 - x[0]) ** power + (x1 - x[1]) ** power + (x2 - x[2]) ** power) / R_power)

        value[:] = - self.normal_component * parabolic


def compute_boundary_geometry_acrn(mesh, ind, facet_domains):
    # Adapted from Womersley.py
    # Convenient information about mesh facets
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    d = mesh.geometry().dim()
    x = SpatialCoordinate(mesh)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0, name="one") * dsi)
    if A == 0:
        return None

    # Compute barycenter by integrating x components over all facets
    c = [assemble(x[i] * dsi) / A for i in range(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = FacetNormal(mesh)
    ni = np.array([assemble(n[i] * dsi) for i in range(d)])
    n_len = np.sqrt(sum([ni[i] ** 2 for i in range(d)]))
    normal = ni / n_len

    # Estimate radius
    r = np.sqrt(A / np.pi)

    return A, c, r, normal
