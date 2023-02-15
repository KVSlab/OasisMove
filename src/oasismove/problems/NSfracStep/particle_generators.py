# flake8: noqa
import numpy as np
from math import pi, sqrt
from itertools import product
from mpi4py import MPI as pyMPI

comm = pyMPI.COMM_WORLD

def circle(center, radius, N, fill=0):
    '''
    With fill ==0:
        put N particles on circle with radius centered at center
    othersise:
        fill circle with smaller concentric circles
    '''
    if fill == 0:
        theta = np.linspace(0, 2.*pi, N, endpoint=False)
        xs = center[0] + radius*np.cos(theta)
        ys = center[1] + radius*np.sin(theta)
        return [np.array([x, y]) for x, y, in zip(xs, ys)]
    else:
        return sum((circle(center, r, N, fill=0)
                    for r in np.linspace(0.1, radius, fill)), [])


class RandomGenerator(object):
    '''
    Fill object by random points.
    '''
    def __init__(self, domain, rule):
        '''
        Domain specifies bounding box for the shape and is used to generate
        points. The rule filter points of inside the bounding box that are
        axctually inside the shape.
        '''
        assert isinstance(domain, list)
        self.domain = domain
        self.rule = rule
        self.dim = len(domain)
        self.rank = comm.Get_rank()

    def generate(self, N, method='full'):
        'Genererate points.'
        assert len(N) == self.dim
        assert method in ['full', 'tensor']
        np.random.seed(2)

        if self.rank == 0:
            # Generate random points for all coordinates
            if method == 'full':
                n_points = np.product(N)
                points = np.random.rand(n_points, self.dim)
                for i, (a, b) in enumerate(self.domain):
                    points[:, i] = a + points[:, i]*(b-a)
            # Create points by tensor product of intervals
            else:
                # Values from [0, 1) used to create points between
                # a, b - boundary
                # points in each of the directiosn
                shifts_i = np.array([np.random.rand(n) for n in N])
                # Create candidates for each directions
                points_i = (a+shifts_i[i]*(b-a)
                            for i, (a, b) in enumerate(self.domain))
                # Cartesian product of directions yield n-d points
                points = (np.array(point) for point in product(*points_i))


            # Use rule to see which points are inside
            points_inside = np.array(list(filter(self.rule, points)))
        else:
            points_inside = None

        points_inside = comm.bcast(points_inside, root=0)

        return points_inside


class RandomRectangle(RandomGenerator):
    def __init__(self, ll, ur):
        # a is lower left, b is upper right
        ax, ay = ll.x(), ll.y()
        bx, by = ur.x(), ur.y()
        assert ax < bx and ay < by
        RandomGenerator.__init__(self, [[ax, bx], [ay, by]], lambda x: True)


class RandomCircle(RandomGenerator):
    def __init__(self, center, radius):
        assert radius > 0
        domain = [[center[0]-radius, center[0]+radius],
                  [center[1]-radius, center[1]+radius]]
        RandomGenerator.__init__(self, domain,
                                 lambda x: sqrt((x[0]-center[0])**2 +
                                                (x[1]-center[1])**2) < radius
                                 )

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from dolfin import Point

    r_rectangle = RandomRectangle(Point(0, 0), Point(1, 5)).generate([100, 100],
                                                            method='tensor')
    r_circle = RandomCircle(Point(0, 0), 1).generate([100, 100])

    for points in [r_rectangle, r_circle]:
        plt.figure()
        plt.scatter(points[:, 0], points[:, 1])
        plt.axis('equal')

    plt.show()
