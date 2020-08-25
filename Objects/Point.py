import math
from Tools.coordinates import *
class Point:
    def __init__(self, **kwargs):
        if ('x' in kwargs) and ('y' in kwargs) and ('z' in kwargs):
            self._x = kwargs['x']
            self._y = kwargs['y']
            self._z = kwargs['z']
        elif ('r' in kwargs) and ('phi' in kwargs) and ('theta' in kwargs):
            self._x, self._y, self._z = spherical_to_cartesian(kwargs['r'],kwargs['phi'],kwargs['theta'])
        elif ('r_cyl' in kwargs) and ('phi' in kwargs) and ('z' in kwargs):
            self._x, self._y, self._z = cylindrical_to_cartesian(kwargs['r_cyl'],kwargs['phi'],kwargs['z'])
        else:
            print(kwargs)
            raise AttributeError('Not enough coordinates to make a point!')

    def x(self):
        return self._x

    def y(self):
        return self._y

    def z(self):
        return self._z

    def r(self):
        return r_from_cartesian(self._x,self._y,self._z)

    def phi(self):
        return phi_from_cartesian(self._x,self._y)

    def theta(self):
        return theta_from_cartesian(self._x,self._y,self._z)

    def r_cyl(self):
        return r_cyl_from_cartesian(self._x,self._y)

    def inverted_point(self):
        x = -1. * self.x()
        y = -1. * self.y()
        z = -1. * self.z()
        return Point(x=x, y=y, z=z)

    def invert(self):
        self._x = -1. * self.x()
        self._y = -1. * self.y()
        self._z = -1. * self.z()

    def __str__(self):
        return 'x:{},y:{},z:{}\nr:{},phi:{},theta:{}'.format(self._x,self._y,self._z,self.r(),self.phi(),self.theta())

def distance(point_a,point_b):
    x = point_a.x() - point_b.x()
    y = point_a.y() - point_b.y()
    z = point_a.z() - point_b.z()
    return math.sqrt(x**2 + y**2 + z**2)
