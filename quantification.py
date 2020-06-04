from Objects.Patterns import SpiralPhyllotaxisPattern, ArchimedeanSpiralUniform, ArchimedeanSpiralNonUniform
from Tools.fibonacci import make_fibonacci
import matplotlib.pyplot as plt

fibonacci = make_fibonacci(13)

averages = []
averages_spir_uni = []
averages_spir_nonuni = []
rsds = []
rsds_spir_uni = []
rsds_spir_nonuni = []
phyllos = []

for fibo in fibonacci:
    print('starting fibo = {}'.format(fibo))
    phyllo = SpiralPhyllotaxisPattern(fibo*34,fibo,alternated_points=False)
    phyllos.append(phyllo)
    averages.append(phyllo.compute_average_distance_between_spokes())
    rsds.append(phyllo.compute_rsd())
    spir_uni = ArchimedeanSpiralUniform(fibo*34,fibo)
    averages_spir_uni.append(spir_uni.compute_average_distance_between_spokes())
    rsds_spir_uni.append(spir_uni.compute_rsd())
    spir_nonuni = ArchimedeanSpiralNonUniform(fibo*34,fibo)
    averages_spir_nonuni.append(spir_nonuni.compute_average_distance_between_spokes())
    rsds_spir_nonuni.append(spir_nonuni.compute_rsd())


plt.scatter(fibonacci,averages,label='phyllotaxis')
plt.scatter(fibonacci,averages_spir_uni,label='uniform Archimedean spiral')
plt.scatter(fibonacci,averages_spir_nonuni,label='nonuniform Archimedean spiral')
plt.xlabel('number of interleaves')
plt.ylabel('average distance to next point')
plt.legend(loc='upper left')
plt.savefig('averages.png')
plt.show()

plt.scatter(fibonacci,rsds,label='phyllotaxis')
plt.scatter(fibonacci,rsds_spir_uni,label='uniform Archimedean spiral')
plt.scatter(fibonacci,rsds_spir_nonuni,label='nonuniform Archimedean spiral')
plt.xlabel('number of interleaves')
plt.ylabel('RSD (%)')
plt.legend(loc='upper left')
plt.savefig('rsds.png')
plt.show()