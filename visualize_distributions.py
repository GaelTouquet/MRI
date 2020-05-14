N_total = 8 * 1597
time_per_acquisition = 0.006

from Objects.Patterns import SpiralPhyllotaxisPattern
article = SpiralPhyllotaxisPattern(1600,55,time_per_acquisition=time_per_acquisition,alternated_points=False)
article.draw(title='Article, basic pattern', alpha=0.8, marker_size=1)

inverted_article = SpiralPhyllotaxisPattern(1600,55,time_per_acquisition=time_per_acquisition,alternated_points=True)
inverted_article.draw(title='Article, alternated basic pattern', alpha=0.8, marker_size=1)

#cardiac phases
from Tools.Heart import Heart
heart = Heart(N_total*time_per_acquisition)
inverted_article.separate_in_time_phases(heart.phases)
inverted_article.draw('Article, alternated basic pattern with cardiac phases', colour='cardiac_phase',alpha=0.8, marker_size=1)

def is_first_cardiac_phase(point):
    return point.cardiac_phase == 0

inverted_article.draw('Article, alternated basic pattern with cardiac phases', colour='cardiac_phase',alpha=0.8, filter_func=is_first_cardiac_phase, marker_size=1)

#now monica's
for n_interleaf in [4,8,16,32]:
    monicas = SpiralPhyllotaxisPattern(N_total,n_interleaf,time_per_acquisition=time_per_acquisition,alternated_points=True)
    monicas.draw(title='N_total = {} , {} interleaves'.format(N_total,n_interleaf), colour='interleaf')

# articles to show n_interleaf does not change pattern, just interleave style
fibo = [1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597]
for n_interleaf in fibo:
    article = SpiralPhyllotaxisPattern(1597,n_interleaf,time_per_acquisition=time_per_acquisition,alternated_points=False)
    article.draw(title='Article N_total = {} , {} interleaf'.format(1597,n_interleaf), colour='interleaf', alpha=0.9, marker_size=1)
