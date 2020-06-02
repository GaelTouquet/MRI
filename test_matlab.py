from Objects.Patterns import SpiralPhyllotaxisPattern

fibo = 1
phyllo = SpiralPhyllotaxisPattern(fibo*34,fibo,time_per_acquisition=0.006,add_readout_ends=False,alternated_points=False)

phyllo.write_to_matlab('./matlab_test.mat','k')