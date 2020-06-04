N_target = 30000
time_per_acquisition = 0.006
fibonacci = [1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711]
from Objects.Patterns import SpiralPhyllotaxisPattern
from Tools.fibonacci import find_best_number_of_points

#cardiac phases
from Tools.Heart import Heart
heart = Heart(N_target*time_per_acquisition * 2,0.05)

#articles to show n_interleaf does not change pattern, just interleave style
def only_first_cardiac_phase(point):
    return point.cardiac_phase == 0

for fibo in fibonacci:
    for add_readout_ends in [True, False]:
        N = find_best_number_of_points(N_target, fibo)
        test = SpiralPhyllotaxisPattern(N*fibo,fibo,time_per_acquisition=time_per_acquisition,alternated_points=True)
        test.separate_in_time_phases(heart.beats)
        test.draw(title='N_total = {} , {} interleaf'.format(N*fibo,fibo), alpha=0.3, marker_size=1, colour='cardiac_phase', filter_func=only_first_cardiac_phase, display=True)