import math

golden_angle = math.pi * (3 - (5**0.5))

# TODO: separate in several modules

### Coordinates

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

### Fibonacci

def make_fibonacci(length):
     if length == 1:
          return [1]
     fibo = [1,2]
     for i in range(length-2):
          fibo.append(fibo[-2] + fibo[-1])
     return fibo

def find_best_number_of_points(Ntarget, fibo):
    i = 0
    previous_number = 0
    new_number = 0
    while (i * fibo) < Ntarget:
        i += 1
        previous_number = new_number
        new_number = i * fibo
    if abs(previous_number - Ntarget) < abs(new_number - Ntarget):
        return i - 1 
    else:
        return i

def find_best_fibo(Ntarget, n_readouts_per_interleaf, fibonacci):
    previous_number = 0
    new_number = 0
    for i in range(len(fibonacci)):
        previous_number = new_number
        new_number = fibonacci[i] * n_readouts_per_interleaf
        if new_number > Ntarget:
            break
    if abs(previous_number - Ntarget) < abs(new_number - Ntarget):
        return fibonacci[i-1]
    else:
        return fibonacci[i]