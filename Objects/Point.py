import math
from Tools.toolbox import cartesian_to_spherical, cartesian_to_cylindrical, spherical_to_cartesian, cylindrical_to_cartesian

class Point:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if hasattr(self, 'x') and hasattr(self, 'y') and hasattr(self, 'z'):
            self.r, self.phi, self.theta = cartesian_to_spherical(self.x, self.y, self.z)
            self.r_cyl, self.phi, self.z = cartesian_to_cylindrical(self.x,self.y,self.z)
        elif hasattr(self, 'r') and hasattr(self, 'phi') and hasattr(self, 'theta'):
            self.x, self.y, self.z = spherical_to_cartesian(self.r,self.phi,self.theta)
            self.r_cyl, self.phi, self.z = cartesian_to_cylindrical(self.x,self.y,self.z)
        elif hasattr(self, 'r_cyl') and hasattr(self, 'phi') and hasattr(self, 'z'):
            self.x, self.y, self.z = cylindrical_to_cartesian(self.r_cyl,self.phi,self.z)
            self.r, self.phi, self.theta = cartesian_to_spherical(self.x, self.y, self.z)
        else:
            print(kwargs)
            raise AttributeError('Not enough coodrinates to make a point!')

    def inverted_point(self):
        x = -1.*self.x
        y = -1*self.y
        z = -1.*self.z
        return Point(x=x, y=y, z=z)

    def __str__(self):
        return 'x:{},y:{},z:{}'.format(self.x,self.y,self.z)

def distance(point_a,point_b):
    x = point_a.x - point_b.x
    y = point_a.y - point_b.y
    z = point_a.z - point_b.z
    return math.sqrt(x**2 + y**2 + z**2)
