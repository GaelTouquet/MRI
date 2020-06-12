import math

golden_angle = math.pi * (3 - (5**0.5))

def x_from_spherical(r,phi,theta):
     return r * math.sin(theta) * math.cos(phi)

def y_from_spherical(r,phi,theta):
     return r * math.sin(theta) * math.sin(phi)

def z_from_spherical(r,theta):
     return r * math.cos(theta)

def spherical_to_cartesian(r,phi,theta):
    return [
         x_from_spherical(r,phi,theta),
         y_from_spherical(r,phi,theta),
         z_from_spherical(r, theta)
    ]

def r_from_cartesian(x,y,z):
     return math.sqrt(x**2 + y**2 + z**2)

def phi_from_cartesian(x,y):
     return math.atan2(y,x)

def theta_from_cartesian(x,y,z):
     return math.atan2(math.sqrt(x**2 + y**2), z)

def cartesian_to_spherical(x, y, z):
     hxy = x**2 + y**2
     return [
          r_from_cartesian(x,y,z),
          phi_from_cartesian(x,y),
          theta_from_cartesian(x,y,z)
     ]

def r_cyl_from_cartesian(x,y):
     return math.sqrt(x**2 + y**2)

def cartesian_to_cylindrical(x, y, z):
     return [
          r_cyl_from_cartesian(x,y),
          phi_from_cartesian(x,y),
          z
     ]

def x_from_cylindrical(r_cyl,phi):
     return r_cyl * math.cos(phi)

def y_from_cylindrical(r_cyl,phi):
     return r_cyl * math.sin(phi)

def cylindrical_to_cartesian(r_cyl, phi, z):
     return [
          x_from_cylindrical(r_cyl,phi),
          y_from_cylindrical(r_cyl,phi),
          z
     ]