import numpy as np
from Objects.Patterns import SpiralPhyllotaxisPattern

fibo = 1
phyllo = SpiralPhyllotaxisPattern(
    8*2584, 2584, n_readouts=272, alternated_points=True)

phyllo.write_to_matlab('./matlab_test.mat','myk',dtype=np.float32)
