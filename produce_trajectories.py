from Objects.Patterns import SpiralPhyllotaxisPattern
from Tools.fibonacci import make_fibonacci, find_best_fibo

fibo = make_fibonacci(15)

n_readouts_per_spoke = 272

for n_spokes_per_interleaf in [8,16]:
    n_interleaves = find_best_fibo(8*272*2594, n_spokes_per_interleaf*272, fibo)
    phyllo = SpiralPhyllotaxisPattern(n_readouts_per_spoke=n_readouts_per_spoke,n_spokes_per_interleaf=n_spokes_per_interleaf,n_interleaves=n_interleaves,rmax=0.5)
    phyllo.write_to_matlab('./matlab_test_full.mat','phyllotaxis_{}'.format(n_spokes_per_interleaf))