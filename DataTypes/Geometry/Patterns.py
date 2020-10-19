import math
import copy
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from MRI.DataTypes.Geometry.Point import Point, distance
from MRI.Tools.fibonacci import golden_angle, golden_angle_2D_1, golden_angle_2D_2

class Pattern:
    """Base class for all patterns. Holds useful methods"""
    def __init__(self):
        self.points = np.array([],dtype=Point)

    def clone(self, filter_func):
        """Create a new pattern which is copied from the pattern but all points pass the given filter function."""
        if filter_func:
            indexes = []
            for i in np.ndindex(self.points.shape):
                if filter_func(self.points[i],i,self.points):
                    indexes.append(i)
            points = np.zeros(len(indexes),dtype=Point)
            for i in range(len(indexes)):
                points[i] = self.points[indexes[i]]
        else:
            points = copy.deepcopy(self.points)
        return CustomPattern(points)

    def get_points(self, filter_func = None, extra_vars = []):
        """Get list of matplotlib-usable points that pass the filter function."""
        points = {}
        points['xs'] = []
        points['ys'] = []
        points['zs'] = []
        for var in extra_vars:
            points[var] = []
        for i in np.ndindex(self.points.shape):
            if (not filter_func) or filter_func(self.points[i]):
                points['xs'].append(self.points[i].x())
                points['ys'].append(self.points[i].y())
                points['zs'].append(self.points[i].z())
                for var in extra_vars:
                    if var == 'first_interleaf':
                        points[var].append('r' if (i[2] == 0) else 'b')
                    else:
                        points[var].append(getattr(self.points[i],var))
        return points

    def separate_in_time_phases(self, time_dict, cycle_name = 'cardiac', phases_names = None, temporal_resolution=None):
        for i in np.ndindex(self.points.shape):
            point = self.points[i]
            phase = 0
            beat = 0
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
        print('selected {} points of {} '.format(len(point_dict['xs']),self.points.size))
        if title:
            plt.title(title)
        if display:
            plt.show()
        if save:
            plt.savefig(save)
        plt.close('all')
        
    def compute_rsd(self):
        for i in np.ndindex(self.points.shape):
            point = self.points[i]
            distances = []
            for j in np.ndindex(self.points.shape):
                other_point = self.points[j]
                if other_point == point:
                    continue
                distances.append(distance(point,other_point))
            point.closest_points_distances = sorted(distances)[:4]
        mu_d = 0
        for i in np.ndindex(self.points.shape):
            for dist in self.points[i].closest_points_distances:
                mu_d += dist
        mu_d = (1/(4*self.n_points)) * mu_d
        rsd = 0
        for i in np.ndindex(self.points.shape):
            for dist in self.points[i].closest_points_distances:
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
        shape = list(self.points.shape)
        shape.append(3)
        savearray = np.zeros(shape,dtype=dtype)
        for i in np.ndindex(self.points.shape):
            x = max(min(self.points[i].x()/norm_factor,0.5),-0.5)
            y = max(min(self.points[i].y()/norm_factor,0.5),-0.5)
            z = max(min(self.points[i].z()/norm_factor,0.5),-0.5)
            savearray[i][:] = x, y, z
        if os.path.isfile(path):
            savedict = loadmat(path)
            savedict[name] = savearray
        else:
            savedict = {name:savearray}
        savemat(path,savedict)

class CustomPattern(Pattern):
    """Pattern with user-defined points."""
    def __init__(self, points):
        self.points = points
        self.n_points = self.points.size

##### Functional patterns

class FunctionalPattern(Pattern):
    """Class for pattern that use a function to derive the position of its points. self.point_function should return a point for each iteration value over n_points."""
    def __init__(self, point_function, n_points, time_per_acquisition=None):
        self.n_points = n_points
        self.point_function = point_function
        self.time_per_acquisition = time_per_acquisition
        if self.time_per_acquisition:
            self.total_time = 0.
        self.update_points()

    def update_points(self):
        self.points = np.zeros(self.n_points,dtype=Point)
        for n in range(self.n_points):
            self.points[n] = self.point_function(n)
            if self.time_per_acquisition:
                self.points[n].t = self.total_time
                self.total_time += self.time_per_acquisition

##### Interleaved patterns

class InterleavedFunctionalPattern(Pattern):
    """Virtual class to hold all the methods useful for the interleaved patterns."""
    def __init__(self, interleaf_function, spoke_function, readout_function, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=None, alternated_spokes=True):
        self.interleaf_function = interleaf_function
        self.spoke_function = spoke_function
        self.readout_function = readout_function
        self.n_readouts_per_spoke = n_readouts_per_spoke
        self.n_spokes_per_interleaf = n_spokes_per_interleaf
        self.n_interleaves = n_interleaves
        self.n_points = self.n_readouts_per_spoke * self.n_spokes_per_interleaf * self.n_interleaves
        self.time_per_acquisition = time_per_acquisition
        self.alternated_spokes = alternated_spokes
        self.update_points()

    def update_points(self):
        self.points = np.zeros((self.n_readouts_per_spoke,self.n_spokes_per_interleaf,self.n_interleaves), dtype=Point)
        if self.time_per_acquisition:
            self.total_time = 0.
        for i_interleaf in range(self.n_interleaves):
            first_spoke = self.interleaf_function(i_interleaf)
            for i_spoke in range(self.n_spokes_per_interleaf):
                first_readout = self.spoke_function(i_interleaf,first_spoke,i_spoke)
                for i_readout in range(self.n_readouts_per_spoke):
                    if self.alternated_spokes and (i_spoke % 2 == 1):
                        k_readout = (self.n_readouts_per_spoke-1) - i_readout
                    else:
                        k_readout = i_readout
                    self.points[i_readout,i_spoke,i_interleaf] = self.readout_function(i_interleaf,first_spoke,i_spoke,first_readout,k_readout)
                    if (self.alternated_spokes and (i_spoke % 2 == 1)) and self.n_readouts_per_spoke==1:
                        self.points[i_readout,i_spoke,i_interleaf].invert()
                    if self.time_per_acquisition and i_readout==0:
                        self.points[i_readout,i_spoke,i_interleaf].t = self.total_time
                        self.total_time += self.time_per_acquisition

    def compute_average_distance_between_spokes(self):
        points_in_first_interleaf = self.points[0,:,0]
        average = 0
        for i in range(len(points_in_first_interleaf)-1):
            average += distance(points_in_first_interleaf[i],points_in_first_interleaf[i+1])
        average = average / (len(points_in_first_interleaf)-1)
        self.average_distance_between_readouts = average
        return average

    def separate_in_time_phases(self, time_dict, cycle_name = 'cardiac', phases_names = None, temporal_resolution=None):
        for i in np.ndindex(self.points.shape):
            point = self.points[i]
            if not hasattr(point,'t'):
                continue
            phase = 0
            beat = 0
            while point.t > time_dict[beat][phase]:
                if phase == len(time_dict[beat])-1:
                    phase = 0
                    beat += 1
                else:
                    phase += 1
            setattr(point, '{}_phase'.format(cycle_name), phase if not phases_names else phases_names[phase])

### Spherical central patterns

class SphericalCentralPattern(InterleavedFunctionalPattern):
    """Class to make and hold spherical patterns."""
    def __init__(self, interleaf_function, spoke_function, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=None, alternated_spokes=False, rmax = 1.):
        self.rmax = rmax
        super().__init__(interleaf_function=interleaf_function, spoke_function=spoke_function, readout_function=self.readout_function, n_readouts_per_spoke=n_readouts_per_spoke, n_spokes_per_interleaf=n_spokes_per_interleaf, n_interleaves=n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes)
    
    def readout_function(self,i_interleaf,first_spoke,i_spoke,first_readout,k_readout):
        if k_readout == 0:
            return first_readout
        else:
            r = first_readout.r()*(1.-2*((k_readout)/(self.n_readouts_per_spoke-1)))
            phi = first_readout.phi()
            theta = first_readout.theta()
            return Point(r=r,phi=phi,theta=theta)

# Specific central patterns

class ArchimedeanSpiralUniform(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.22898 and https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.20128"""
    def __init__(self, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=None, alternated_spokes=False, rmax = 1.):
        super().__init__(interleaf_function=self.interleaf_function, spoke_function=self.spoke_function, n_readouts_per_spoke=n_readouts_per_spoke, n_spokes_per_interleaf=n_spokes_per_interleaf, n_interleaves=n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes, rmax = rmax)
    
    def spoke_function(self, i_interleaf,first_spoke,i_spoke):
        n = i_spoke*self.n_interleaves + i_interleaf
        z = self.rmax * (1 - (n/(self.n_interleaves*self.n_spokes_per_interleaf)))
        phi = (math.sqrt(2*self.n_interleaves*self.n_spokes_per_interleaf*math.pi)*math.asin(z))
        x = math.cos(phi)*math.sqrt(1-(z*z))
        y = math.sin(phi)*math.sqrt(1-(z*z))
        return Point(x=x,y=y,z=z)

    def interleaf_function(self, i_interleaf):
        return None

class ArchimedeanSpiralNonUniform(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.22898 and https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.20128"""
    def __init__(self, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=None, alternated_spokes=False, rmax = 1.):
        super().__init__(interleaf_function=self.interleaf_function, spoke_function=self.spoke_function, n_readouts_per_spoke=n_readouts_per_spoke, n_spokes_per_interleaf=n_spokes_per_interleaf, n_interleaves=n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes, rmax = rmax)
    
    def spoke_function(self, i_interleaf,first_spoke,i_spoke):
        if i_interleaf == 0:
            z = self.rmax * (1 - (i_spoke/self.n_spokes_per_interleaf))
            phi = (math.sqrt((2*self.n_spokes_per_interleaf*math.pi)/self.n_interleaves)*math.asin(z)) + ((2*i_interleaf*math.pi)/self.n_interleaves)
            x = math.cos(phi)*math.sqrt((self.rmax*self.rmax)-(z*z))
            y = math.sin(phi)*math.sqrt((self.rmax*self.rmax)-(z*z))
            return Point(x=x,y=y,z=z)
        else:
            r = self.points[0,i_spoke,0].r()
            phi = self.points[0,i_spoke,0].phi() + (i_interleaf*2*math.pi/self.n_interleaves)
            theta = self.points[0,i_spoke,0].theta()
            return Point(r=r,phi=phi,theta=theta)

    def interleaf_function(self, i_interleaf):
        return None

class SpiralPhyllotaxisPattern(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/10.1002/mrm.22898"""
    def __init__(self, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=None, alternated_spokes=True, rmax=1.0):
        super().__init__(interleaf_function=self.interleaf_function, spoke_function=self.spoke_function, n_readouts_per_spoke=n_readouts_per_spoke, n_spokes_per_interleaf=n_spokes_per_interleaf, n_interleaves=n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes, rmax=rmax)

    def spoke_function(self, i_interleaf,first_spoke,i_spoke):
        n = i_spoke*self.n_interleaves + i_interleaf
        r = self.rmax
        theta = (math.pi/2) * math.sqrt(n/(self.n_interleaves*self.n_spokes_per_interleaf))
        phi = n * golden_angle
        return Point(r=r, theta=theta, phi=phi)

    def interleaf_function(self, i_interleaf):
        return None

class MultiGoldenRadial(SphericalCentralPattern):
    """https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.21837"""
    def __init__(self, n_readouts_per_spoke, n_spokes, time_per_acquisition=None, alternated_spokes=True, rmax=1.0):
        super().__init__(interleaf_function=lambda x: None, spoke_function=self.spoke_function, n_readouts_per_spoke=n_readouts_per_spoke, n_spokes_per_interleaf=n_spokes, n_interleaves=1, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes, rmax=rmax)

    def spoke_function(self, i_interleaf,first_spoke,i_spoke):
        z = (i_spoke * golden_angle_2D_1) % 1
        phi = i_spoke * golden_angle_2D_2
        r = self.rmax
        r_cyl = math.sqrt(r*r - z*z)
        return Point(r_cyl=r_cyl, phi=phi, z=z)

# Interleaved pattern from matlab format

class InterleavedMatLabPattern(InterleavedFunctionalPattern):
    """Class to load a pattern from a saved MatLab file."""
    def __init__(self, path, time_per_acquisition=None, alternated_spokes=True):
        matlab_contents = loadmat(path)
        matlab_pattern = matlab_contents['k']
        del(matlab_contents)
        shape = matlab_pattern.shape
        n_readouts_per_spoke = shape[0]
        n_spokes_per_interleaf = shape[1]
        n_interleaves = shape[2]
        self.points = np.zeros((n_readouts_per_spoke,n_spokes_per_interleaf,n_interleaves),dtype=Point)
        for i_readout in range(n_readouts_per_spoke):
            for i_spoke in range(n_spokes_per_interleaf):
                for i_interleaf in range(n_interleaves):
                    point_coordinates = matlab_pattern[i_readout][i_spoke][i_interleaf]
                    self.points[i_readout,i_spoke,i_interleaf] = Point(x=point_coordinates[0],y=point_coordinates[1],z=point_coordinates[2])
        del(matlab_pattern)
        super().__init__(None, None, None, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes)

    def update_points(self):
        if self.time_per_acquisition:
            self.total_time = 0.
            for i_interleaf in range(self.n_interleaves):
                for i_spoke in range(self.n_spokes_per_interleaf):
                    for i_readout in range(self.n_readouts_per_spoke):
                        self.points[i_readout,i_spoke,i_interleaf].t = self.total_time + i_readout*(self.time_per_acquisition/self.n_readouts_per_spoke)
                    self.total_time += self.time_per_acquisition

### 2D Stacks Patterns

class StackPattern(InterleavedFunctionalPattern):
    """Base class for 2D Stacked patterns."""
    def __init__(self, spoke_function, readout_function, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, total_height=1., time_per_acquisition=None, alternated_spokes=True):
        self.total_height = total_height
        super().__init__(self.interleaf_function, spoke_function, readout_function, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes)

    def interleaf_function(self, i_interleaf):
        return Point(x=0,y=0,z=self.total_height*(1-2*((i_interleaf)/(self.n_interleaves-1))))

# Specific stack patterns

class StackOfStarsPattern(StackPattern):
    """Stack of stars pattern."""
    def __init__(self, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, total_height=1., rmax=1., time_per_acquisition=None, alternated_spokes=True):
        self.rmax = rmax
        super().__init__(self.spoke_function, self.readout_function, n_readouts_per_spoke, n_spokes_per_interleaf, n_interleaves, total_height=total_height, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes)

    def spoke_function(self, i_interleaf,first_spoke,i_spoke):
        z = first_spoke.z()
        r_cyl = self.rmax
        phi = i_spoke*math.pi/self.n_spokes_per_interleaf
        return Point(r_cyl=r_cyl,phi=phi,z=z)

    def readout_function(self,i_interleaf,first_spoke,i_spoke,first_readout,k_readout):
        if k_readout==0:
            return first_readout
        else:
            z = first_readout.z()
            phi = first_readout.phi()
            r_cyl = first_readout.r_cyl()*(1.-2*((k_readout)/(self.n_readouts_per_spoke-1)))
            return Point(r_cyl=r_cyl,phi=phi,z=z)