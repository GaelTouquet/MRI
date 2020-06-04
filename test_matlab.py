from Objects.Patterns import SpiralPhyllotaxisPattern

fibo = 1
phyllo = SpiralPhyllotaxisPattern(
    8*2584, 2584, n_readouts=100, time_per_acquisition=0.006, alternated_points=True)

# phyllo.write_to_matlab('./matlab_test.mat','k')
