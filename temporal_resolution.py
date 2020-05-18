from Objects.Patterns import SpiralPhyllotaxisPattern
from Tools.toolbox import make_fibonacci, find_best_fibo
from Tools.Heart import Heart

temporal_resolutions = [0.05,0.06,0.12]
heart_beat_lengths = [0.5,1,1.5]
N_readouts_target = 30000
time_per_readout = 0.006
numbers_of_projections_per_interleaf = [4,8,12,16,32,64,1000,5000,15000]

fibonacci = make_fibonacci(100)

beating_time = 2 * N_readouts_target * time_per_readout # to make sure we cover total readout time
hearts = []
for heart_beat_length in heart_beat_lengths:
    hearts.append(Heart(beating_time,n_phases=10,average_bps=heart_beat_length,bps_rms=0.05,bps_lower_limit=heart_beat_length-0.05,bps_higher_limit=heart_beat_length+0.05))

def only_first_cardiac_phase(point):
    return point.cardiac_phase == 0

for temporal_resolution in temporal_resolutions:
    for heart in hearts:
        for number_of_projections_per_interleaf in numbers_of_projections_per_interleaf:
            print('tres {}; heartbeat {}; n proj {}'.format(temporal_resolution,heart.average_bps,number_of_projections_per_interleaf))
            n_interleaf = find_best_fibo(N_readouts_target,number_of_projections_per_interleaf,fibonacci)
            phyllo = SpiralPhyllotaxisPattern(n_interleaf*number_of_projections_per_interleaf,n_interleaf,time_per_acquisition=time_per_readout,add_readout_ends=True,alternated_points=True)
            phyllo.separate_in_time_phases(heart.phases, temporal_resolution=temporal_resolution)
            phyllo.draw(title='N interleaves = {} ; {} readouts per interleaf'.format(phyllo.n_interleaf,phyllo.n_points/phyllo.n_interleaf),alpha=0.3, marker_size=1, colour='cardiac_phase', filter_func=only_first_cardiac_phase, save='pres_18_05_2020/tres{}_beat{}_proj{}.png'.format(temporal_resolution,heart.average_bps,number_of_projections_per_interleaf))