import math
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from Objects.Point import Point, distance
from Tools.coordinates import golden_angle

class Pattern:
    """Base class for all patterns."""
    def __init__(self):
        self.points = np.array([],dtype=Point)

    def clone(self, filter_func):
        """Create a new pattern which is copied from the pattern but all points pass the given filter function."""
        new_pattern = copy.deepcopy(self)
        new_pattern.points = np.array([],dtype=Point)
        for point in self.points:
            if filter_func(point):
                new_pattern.points = np.append(new_pattern.points,point)

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
                points['xs'].append(point.x())
                points['ys'].append(point.y())
                points['zs'].append(point.z())
                for var in extra_vars:
                    if var == 'first_interleaf':
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

    def draw(self, title = None, filter_func=None, colour = None, marker_size=0.1,alpha=0.6, display=False, save=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if colour:
            point_dict = self.get_points(extra_vars=[colour], filter_func=filter_func)
            colour = point_dict[colour]
        else:
            if hasattr(self.points[0], 'interleaf'):
                point_dict = self.get_points(filter_func=filter_func, extra_vars=['first_interleaf'])
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
            for other_point in self.points: #TODO use np.vectorize(distance) as the function to use here? test if gain in time?
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

    def write_to_matlab(self,path,name,dtype=np.float64):
        norm_factor = 2
        if hasattr(self,'rmax'):
            norm_factor = 2*self.rmax
        savearray = np.array([[point.x()/norm_factor,point.y()/norm_factor,point.z()/norm_factor] for point in self.points],dtype=dtype)
        savedict = {name:savearray}
        savemat(path,savedict)

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
        self.points = np.zeros(self.n_points,dtype=Point)
        for n in range(self.n_points):
            self.points[n] = self.point_function(n)

class InterleavedFunctionalPattern(FunctionalPattern):
    """Virtual class to hold all the methods useful for the interleaved patterns."""
    #TODO make interleav init in __init__ of this class
    def compute_average_distance_between_spokes(self):
        points_in_first_interleaf = self.interleaves[0]
        average = 0
        for i in range(len(points_in_first_interleaf)-1):
            average += distance(points_in_first_interleaf[i][0],points_in_first_interleaf[i+1][0])
        average = average / (len(points_in_first_interleaf)-1)
        self.average_distance_between_readouts = average
        return average

class SphericalCentralPattern(InterleavedFunctionalPattern):
    """Class to make and hold spherical patterns."""
    def __init__(self, point_function, n_points, time_per_acquisition=None, n_readouts=1, rmax = 1.):
        self.n_readouts = n_readouts
        self.rmax = rmax
        super().__init__(point_function, n_points, time_per_acquisition=time_per_acquisition)

    def update_points(self):
        self.points = np.empty((self.n_points*self.n_readouts),dtype=Point)
        for n in range(self.n_points):
            point = self.point_function(n)
            self.points[n*self.n_readouts] = point
            if self.n_readouts:
                for i in range(self.n_readouts-1):
                    new_point = Point(r=point.r()*(1.-2*((i+1)/(self.n_readouts-1))),phi=point.phi(),theta=point.theta())
                    self.points[(n*self.n_readouts) + (i+1)] = new_point
                    if hasattr(self,'interleaves'):
                        new_point.interleaf = point.interleaf
                        self.interleaves[point.interleaf][-1].append(new_point)

class ArchimedeanSpiralUniform(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.22898 and https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.20128"""
    def __init__(self, n_points, n_interleaves, time_per_acquisition=None, n_readouts=1, rmax=1.0):
        self.n_interleaf = n_interleaves
        self.n_points_per_interleaf = n_points // self.n_interleaf
        self.interleaves = {}
        for k in range(self.n_interleaf):
            self.interleaves[k] = []
        super().__init__(self.point_function, n_points, time_per_acquisition=time_per_acquisition, n_readouts=n_readouts, rmax=rmax)
    
    def point_function(self,n):
        # p = n % self.n_points_per_interleaf
        # i = n // self.n_points_per_interleaf
        z = self.rmax * (1 - (n/self.n_points))
        # z = (2*p - (self.n_points_per_interleaf + 1))/self.n_points_per_interleaf
        phi = (math.sqrt(2*self.n_points*math.pi)*math.asin(z))
        x = math.cos(phi)*math.sqrt(1-(z*z))
        y = math.sin(phi)*math.sqrt(1-(z*z))
        point = Point(x=x,y=y,z=z)
        m = n % self.n_interleaf
        self.interleaves[m].append([point])
        point.interleaf = m
        return point

class ArchimedeanSpiralNonUniform(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.22898 and https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.20128"""
    def __init__(self, n_points, n_interleaves, time_per_acquisition=None, n_readouts=1, rmax=1.0):
        self.n_interleaf = n_interleaves
        self.n_points_per_interleaf = n_points // self.n_interleaf
        self.interleaves = {}
        for k in range(self.n_interleaf):
            self.interleaves[k] = []
        super().__init__(self.point_function, n_points, time_per_acquisition=time_per_acquisition, n_readouts=n_readouts, rmax=rmax)
    
    def point_function(self,n):
        p = n % self.n_points_per_interleaf
        i = n // self.n_points_per_interleaf
        if i == 0:
            z = self.rmax * (1 - (p/self.n_points_per_interleaf))
            # z = (p - (self.n_interleaf + 0.5))/self.n_points_per_interleaf
            phi = (math.sqrt((2*self.n_points_per_interleaf*math.pi)/self.n_interleaf)*math.asin(z)) + ((2*i*math.pi)/self.n_interleaf)
            x = math.cos(phi)*math.sqrt((self.rmax*self.rmax)-(z*z))
            y = math.sin(phi)*math.sqrt((self.rmax*self.rmax)-(z*z))
            point = Point(x=x,y=y,z=z)
        else:
            r = self.interleaves[0][p][0].r()
            phi = self.interleaves[0][p][0].phi() + (i*2*math.pi/self.n_interleaf)
            theta = self.interleaves[0][p][0].theta ()
            point = Point(r=r,phi=phi,theta=theta)
        self.interleaves[i].append([point])
        point.interleaf = i
        return point

class SpiralPhyllotaxisPattern(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/10.1002/mrm.22898"""
    def __init__(self, n_points, n_interleaf, time_per_acquisition=None, n_readouts=1, alternated_points=True, rmax = 1.):
        self.n_interleaf = n_interleaf
        self.interleaves = {}
        self.alternated_points = alternated_points
        for k in range(self.n_interleaf):
            self.interleaves[k] = []
        super().__init__(self.point_function, n_points, time_per_acquisition=time_per_acquisition, n_readouts=n_readouts, rmax=rmax)
        if self.alternated_points:
            self.alternate_points()

    def alternate_points(self):
        for i, interleaf in self.interleaves.items():
            for k in range(len(interleaf)):
                if k % 2 == 1:
                    interleaf[k].reverse()

    def point_function(self, n):
        r = self.rmax
        theta = (math.pi/2) * math.sqrt(n/self.n_points)
        phi = n * golden_angle
        new_point = Point(r=r, theta=theta, phi=phi)
        k = n % self.n_interleaf
        self.interleaves[k].append([new_point])
        new_point.interleaf = k
        return new_point

    def time_points(self):
        for k_interleaf in self.interleaves:
            interleaf = self.interleaves[k_interleaf]
            for spoke in interleaf:
                for point in spoke:
                    point.t = self.total_time
                    self.total_time += self.time_per_acquisition
                

class MatLabPattern(Pattern):
    """Class to load a pattern from a saved MatLab file."""
    def __init__(self, path, collection_name):
        super().__init__()
        self.load_from_matlab(path, collection_name)

    def load_from_matlab(self, path, collection_name):
        matlab_contents = loadmat(path)
        matlab_pattern = matlab_contents[collection_name]
        del(matlab_contents)
        self.points = []
        shape = matlab_pattern.shape
        for i_readout in range(shape[0]):
            for i_spoke in range(shape[1]):
                for i_interleaf in range(shape[2]):
                    point_coordinates = matlab_pattern[i_readout][i_spoke][i_interleaf]
                    self.points.append(Point(x=point_coordinates[0],y=point_coordinates[1],z=point_coordinates[2]))
        del(matlab_pattern)