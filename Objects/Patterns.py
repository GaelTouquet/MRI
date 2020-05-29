import math
import copy
import matplotlib.pyplot as plt
from Objects.Point import Point, distance, distance
from Tools.coordinates import golden_angle

class Pattern:
    """Base class for all patterns."""
    def __init__(self):
        self.points = []

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

    def separate_in_time_phases(self, time_dict, cycle_name = 'cardiac', phases_names = None, temporal_resolution=None):
        phase = 0
        beat = 0
        for point in sorted(self.points, key = lambda point: point.t):
            while point.t > time_dict[beat][phase]:
                if phase == len(time_dict[beat])-1:
                    phase = 0
                    beat += 1
                else:
                    phase += 1
            setattr(point, '{}_phase'.format(cycle_name), phase if not phases_names else phases_names[phase])

    def draw(self, title = None, filter_func=None, colour = 'first_interleave', marker_size=0.1,alpha=0.6, display=False, save=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if colour:
            point_dict = self.get_points(extra_vars=[colour], filter_func=filter_func)
            colour = point_dict[colour]
        else:
            point_dict = self.get_points(filter_func=filter_func)
        s = [marker_size for x in point_dict['xs']]
        ax.scatter(point_dict['xs'],point_dict['ys'],point_dict['zs'], s=s,c=colour,marker='.',alpha=alpha)
        print('selected {} points of {} '.format(len(point_dict['xs']),len(self.points)))
        if title:
            plt.title(title)
        if display:
            plt.show()
        if save:
            plt.savefig(save)
        plt.close('all')
        
    def compute_rsd(self):
        for point in self.points:
            distances = []
            for other_point in self.points:
                if other_point == point:
                    continue
                distances.append(distance(point,other_point))
            point.closest_points_distances = sorted(distances)[:4]
        mu_d = 0
        for point in self.points:
            for dist in point.closest_points_distances:
                mu_d += dist
        mu_d = (1/(4*self.n_points)) * mu_d
        rsd = 0
        for point in self.points:
            for dist in point.closest_points_distances:
                rsd += ((dist - mu_d) * (dist - mu_d))
        rsd = (1/(4*self.n_points)) * rsd
        rsd = math.sqrt(rsd)
        rsd = (100/mu_d) * rsd
        self.rsd = rsd
        return rsd

class CustomPattern(Pattern):
    """Pattern with user-defined points."""
    def __init__(self, points):
        self.points = points

class FunctionalPattern(Pattern):
    """Base class for pattern that use a function to derive the position of its points. self.point_function should return a point for each iteration value over n_points."""
    def __init__(self, point_function, n_points, time_per_acquisition=None):
        self.n_points = n_points
        self.point_function = point_function
        super().__init__()
        self.update_points()
        self.time_per_acquisition = time_per_acquisition
        if self.time_per_acquisition:
            self.total_time = 0
            self.time_points()

    def update_points(self):
        self.points = []
        for n in range(self.point_function):
            point = self.point_function(n)
            self.points.append(point)

    def time_points(self):
        for i, interleaf in self.interleaves.items():
            for k in range(len(interleaf)):
                interleaf[k].t = self.total_time
                if (not hasattr(self, 'add_readout_ends')) or (not self.add_readout_ends) or (self.add_readout_ends and k % 2 == 1):
                    self.total_time += self.time_per_acquisition

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
                if hasattr(point,'interleaf'):
                    inverted_point.interleaf = point.interleaf
                    self.interleaves[point.interleaf].insert(self.interleaves[point.interleaf].index(point)+1,inverted_point)
            self.points.append(point)
            if self.add_readout_ends:
                self.points.append(inverted_point)

class SpiralPhyllotaxisPattern(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/10.1002/mrm.22898"""
    def __init__(self, n_points, n_interleaf, time_per_acquisition=None, add_readout_ends=False, alternated_points=True, rmax = 1.):
        self.n_interleaf = n_interleaf
        self.interleaves = {}
        self.alternated_points = alternated_points
        for k in range(self.n_interleaf):
            self.interleaves[k] = []
        super().__init__(self.point_function, n_points, time_per_acquisition=time_per_acquisition, add_readout_ends=add_readout_ends, rmax=rmax)
        if self.alternated_points:
            self.alternate_points()

    def alternate_points(self):
        for i, interleaf in self.interleaves.items():
            for k in range(len(interleaf)):
                if self.add_readout_ends:
                    if k % 4 in [2,3]: #swap readout start and readout end for every second pair
                        interleaf[k].invert()
                else:
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

    def compute_average_distance_between_readouts(self):
        points_in_first_interleaf = self.interleaves[0]
        average = 0
        for i in range(len(points_in_first_interleaf)-1):
            average += distance(points_in_first_interleaf[i],points_in_first_interleaf[i+1])
        average = average / (len(points_in_first_interleaf)-1)
        self.average_distance_between_readouts = average
        return average

class Golden3DRadialPatern(SphericalCentralPattern):
    def __init__(self, point_function, n_points, time_per_acquisition=None, add_readout_ends=False, rmax=1.0, golden_1=0.4656, golden_2=0.6823):
        self.golden_1 = golden_1
        self.golden_2 = golden_2
        super().__init__(point_function, n_points, time_per_acquisition=time_per_acquisition, add_readout_ends=add_readout_ends, rmax=rmax)

    def point_function(self, n):
        r = self.rmax
        phi = n * self.golden_1
        z = ( (n * self.golden_2) % 1 ) * self.rmax
        r_cyl = r * math.sin()
        new_point = Point()
