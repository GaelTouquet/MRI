import random

class Heart:
    """Class that computes and holds the information of a heart beating."""
    def __init__(self, total_beating_time, n_phases = 20, average_bps = 1., bps_rms = 0.25, bps_lower_limit = 0.5, bps_higher_limit = 1.5, starting_time = -0.42):
        self.n_phases = n_phases
        self.average_bps = average_bps
        self.bps_rms = bps_rms
        self.bps_lower_limit = bps_lower_limit
        self.bps_higher_limit = bps_higher_limit
        self.total_beating_time = total_beating_time
        self.starting_time = starting_time
        self.play_your_funky_music()

    def beat(self):
        beat_time = random.gauss(self.average_bps,self.bps_rms)
        while beat_time < self.bps_lower_limit or beat_time > self.bps_higher_limit:
            beat_time = random.gauss(self.average_bps,self.bps_rms)
        return beat_time

    def play_your_funky_music(self):
        self.beats = {}
        self.phases = {}
        for n in range(self.n_phases):
            self.phases[n] = []
        self.n_beats_so_far = 0
        self.time_so_far = self.starting_time
        while (self.time_so_far<self.total_beating_time):
            beat_time = self.beat()
            self.beats[self.n_beats_so_far] = beat_time
            self.make_phases_from_beat(beat_time)
            self.time_so_far += beat_time
            self.n_beats_so_far += 1
        
    def make_phases_from_beat(self,beat_time):
        phase_length = beat_time/self.n_phases
        for n in range(self.n_phases):
            self.phases[n].append(self.time_so_far+(n*phase_length))
