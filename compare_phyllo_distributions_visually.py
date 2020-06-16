from Objects.Patterns import SpiralPhyllotaxisPattern, ArchimedeanSpiralUniform, ArchimedeanSpiralNonUniform
from Tools.Heart import Heart

n_spokes = 5
n_interleaves = 8

phyllo = SpiralPhyllotaxisPattern(n_readouts_per_spoke=1,n_spokes_per_interleaf=n_spokes,n_interleaves=n_interleaves,alternated_spokes=False)
phyllo.draw(colour='first_interleaf',marker_size=10,display=True)

arch_uni = ArchimedeanSpiralUniform(n_readouts_per_spoke=1,n_spokes_per_interleaf=n_spokes,n_interleaves=n_interleaves)
arch_uni.draw(colour='first_interleaf',marker_size=10,display=True)

arch_nonuni = ArchimedeanSpiralNonUniform(n_readouts_per_spoke=1,n_spokes_per_interleaf=n_spokes,n_interleaves=n_interleaves)
arch_nonuni.draw(colour='first_interleaf',marker_size=10,display=True)

### temporal resolution
temporal_resolution = 0.06
beating_time = 2 * n_spokes*n_interleaves * temporal_resolution

heart = Heart(beating_time, temporal_phase_resolution=temporal_resolution, average_bps=1., bps_rms=0.0)

phyllo = SpiralPhyllotaxisPattern(n_readouts_per_spoke=1,n_spokes_per_interleaf=n_spokes,n_interleaves=n_interleaves,alternated_spokes=True,time_per_acquisition=0.006)
phyllo.separate_in_time_phases(heart.beats)
phyllo.draw(title='N interleaves = {} ; {} readouts per interleaf ; Ntot = {}'.format(phyllo.n_interleaves, phyllo.n_points/phyllo.n_interleaves, phyllo.n_points), alpha=0.3, marker_size=10, colour='cardiac_phase', display=True)