import math
import copy
import matplotlib.pyplot as plt
from Objects.Point import Point, distance
from Tools.toolbox import golden_angle

class Pattern:
    """Base class for all patterns."""
    def __init__(self, time_per_acquisition=None):
        self.points = []
        if time_per_acquisition:
            self.time_per_acquisition = time_per_acquisition
            self.total_time = 0

    def clone(self, filter_func):
        """Create a new pattern which is copied from the pattern but all points pass the given filter function."""
        new_pattern = copy.deepcopy(self)
        new_pattern.points = []
        for point in self.points:
            if filter_func(point):
                new_pattern.points.append(point)

    def get_points(self, filter_func = None, extra_vars = []):
        """Get list of matplotlib-usable points that pass the filter function."""
        points = {}
        points['xs'] = []
        points['ys'] = []
        points['zs'] = []
        for var in extra_vars:
            points[var] = []
        for point in self.points:
            if (not filter_func) or filter_func(point):
                points['xs'].append(point.x)
                points['ys'].append(point.y)
                points['zs'].append(point.z)
                for var in extra_vars:
                    if var == 'first_interleave':
                        points[var].append('r' if point.interleaf == 0 else 'b')
                    else:
                        points[var].append(getattr(point,var))

        return points

    def separate_in_time_phases(self, time_dict, cycle_name = 'cardiac', phases_names = None):
        n_phase = len(time_dict)
        phase = 0
        beat = 0
        for point in self.points:
            while point.t > time_dict[phase][beat]:
                if phase == n_phase - 1:
                    phase = 0
                    beat += 1
                else:
                    phase += 1
            setattr(point, '{}_phase'.format(cycle_name), phase if not phases_names else phases_names[phase])


    def draw(self, title = None, filter_func=None, colour = 'first_interleave', marker_size=0.1,alpha=0.6):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if colour:
            point_dict = self.get_points(extra_vars=[colour], filter_func=filter_func)
            colour = point_dict[colour]
        else:
            point_dict = self.get_points(filter_func=filter_func)
        s = [marker_size for x in point_dict['xs']]
        ax.scatter(point_dict['xs'],point_dict['ys'],point_dict['zs'], s=s,c=colour,marker='.',alpha=alpha)
        if title:
            plt.title(title)
        plt.show()

class CustomPattern(Pattern):
    """Pattern with user-defined points."""
    def __init__(self, points):
        self.points = points

class FunctionalPattern(Pattern):
    """Base class for pattern that use a function to derive the position of its points. self.point_function should return a point for each iteration value over n_points."""
    def __init__(self, point_function, n_points, time_per_acquisition=None):
        self.n_points = n_points
        self.point_function = point_function
        super().__init__(time_per_acquisition=time_per_acquisition)
        self.update_points()

    def update_points(self):
        self.points = []
        for n in range(self.point_function):
            point = self.point_function(n)
            if hasattr(self, 'time_per_acquisition'):
                point.t = self.total_time
                self.total_time += self.time_per_acquisition
            self.points.append(point)

class SphericalCentralPattern(FunctionalPattern):
    """Class to make and hold spherical patterns."""
    def __init__(self, point_function, n_points, time_per_acquisition=None, add_readout_ends=False, rmax = 1.):
        self.add_readout_ends = add_readout_ends
        self.rmax = rmax
        super().__init__(point_function, n_points, time_per_acquisition=time_per_acquisition)

    def update_points(self):
        self.points = []
        for n in range(self.n_points):
            point = self.point_function(n)
            if self.add_readout_ends:
                inverted_point = point.inverted_point()
            if hasattr(self, 'time_per_acquisition'):
                point.t = self.total_time
                if self.add_readout_ends:
                    inverted_point.t = self.total_time
                    if hasattr(point,'interleaf'):
                        inverted_point.interleaf = point.interleaf
                self.total_time += self.time_per_acquisition
            self.points.append(point)
            if self.add_readout_ends:
                self.points.append(inverted_point)

class SpiralPhyllotaxisPattern(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/10.1002/mrm.22898"""
    def __init__(self, n_points, n_interleaf, time_per_acquisition=None, add_readout_ends=False, alternated_points=True, rmax = 1.):
        self.n_interleaf = n_interleaf
        self.interleaves = {}
        for k in range(self.n_interleaf):
            self.interleaves[k] = []
        super().__init__(self.point_function, n_points, time_per_acquisition=time_per_acquisition, add_readout_ends=add_readout_ends, rmax=rmax)
        if alternated_points:
            self.alternate_points()

    def alternate_points(self):
        for i, interleaf in self.interleaves.items():
            for k in range(len(interleaf)):
                if k % 2 == 1:
                    interleaf[k].invert()

    def point_function(self, n):
        r = self.rmax
        theta = (math.pi/2) * math.sqrt(n/self.n_points)
        phi = n * golden_angle
        new_point = Point(r=r, theta=theta, phi=phi)
        k = n % self.n_interleaf
        self.interleaves[k].append(new_point)
        new_point.interleaf = k
        return new_point