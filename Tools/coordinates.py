import math

golden_angle = math.pi * (3 - (5**0.5))

def spherical_to_cartesian(r,phi,theta):
    return [
         r * math.sin(theta) * math.cos(phi),
         r * math.sin(theta) * math.sin(phi),
         r * math.cos(theta)
    ]

def cartesian_to_spherical(x, y, z):
     hxy = x**2 + y**2
     return [
          math.sqrt(hxy + z**2),
          math.atan2(y,x),
          math.atan2(z, math.sqrt(hxy**2))
     ]

def cartesian_to_cylindrical(x, y, z):
     return [
          math.sqrt(x**2 + y**2),
          math.atan2(y,x),
          z
     ]

def cylindrical_to_cartesian(r_cyl, phi, z):
     return [
          r_cyl * math.cos(phi),
          r_cyl * math.sin(phi),
          z
     ]
