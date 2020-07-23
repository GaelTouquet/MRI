from Objects.Patterns import SpiralPhyllotaxisPattern
from Tools.fibonacci import find_best_number_of_points

n_readouts_per_spoke = 272

for n_spokes_per_interleaf in [8,16]:
    n_interleaves = find_best_number_of_points(8*272*2594)
    phyllo = SpiralPhyllotaxisPattern(n_readouts_per_spoke=n_readouts_per_spoke,n_spokes_per_interleaf=n_spokes_per_interleaf,n_interleaves=n_interleaves,rmax=0.5)
    phyllo.write_to_matlab('./matlab_test_full.mat','phyllotaxis_{}'.format(n_spokes_per_interleaf))