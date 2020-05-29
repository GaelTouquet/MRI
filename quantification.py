from Objects.Patterns import SpiralPhyllotaxisPattern
from Tools.fibonacci import make_fibonacci
import matplotlib.pyplot as plt

fibonacci = make_fibonacci(8)

averages = []
rsds = []
phyllos = []

for fibo in fibonacci:
    print('starting fibo = {}'.format(fibo))
    phyllo = SpiralPhyllotaxisPattern(fibo*34,fibo,time_per_acquisition=0.006,add_readout_ends=False,alternated_points=False)
    phyllos.append(phyllo)
    averages.append(phyllo.compute_average_distance_between_readouts())
    rsds.append(phyllo.compute_rsd())

plt.scatter(fibonacci,averages)
plt.xlabel('number of interleaves')
plt.ylabel('average distance to next point')
plt.savefig('19_09_2020/averages.png')
plt.show()

plt.scatter(fibonacci,rsds)
plt.xlabel('number of interleaves')
plt.ylabel('RSD (%)')
plt.savefig('19_09_2020/rsds.png')
plt.show()

for phyllo in phyllos:
    phyllo.distances_distribution = []
    for point in phyllo.points:
        for distance in point.closest_points_distances:
            phyllo.distances_distribution.append(distance)
