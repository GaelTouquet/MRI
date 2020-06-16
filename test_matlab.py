from Objects.Patterns import SpiralPhyllotaxisPattern, InterleavedMatLabPattern

phyllo = SpiralPhyllotaxisPattern(n_readouts_per_spoke=272,n_spokes_per_interleaf=8,n_interleaves=2584,alternated_spokes=True,rmax=0.5)
phyllo.write_to_matlab('./matlab_test.mat','myk')

# matlab = InterleavedMatLabPattern('ph_11072019_2035184_5_2_wip_phyllo_classicV4_k.mat','k',time_per_acquisition=0.06)
# matlab.draw(filter_func=is_on_sphere,display=True)

#differences: mainly on phi, start at different positions, rotating opposite way ; also alteranting evens for me, odds for monica