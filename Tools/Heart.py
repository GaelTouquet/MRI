import random

class Heart:
    """Class that computes and holds the information of a heart beating."""
    def __init__(self, total_beating_time, temporal_phase_resolution, average_bps = 1., bps_rms = 0.25, bps_lower_limit = 0.5, bps_higher_limit = 1.5, starting_time = -0.42):
        self.temporal_phase_resolution = temporal_phase_resolution
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
        random.seed(0xdeadbeef)
        self.beats = {}
        self.n_beats_so_far = 0
        self.time_so_far = self.starting_time
        while (self.time_so_far<self.total_beating_time):
            beat_time = self.beat()
            self.make_phases_from_beat(beat_time)
            self.time_so_far += beat_time
            self.n_beats_so_far += 1
        
    def make_phases_from_beat(self,beat_time):
        time_within_beat = 0
        self.beats[self.n_beats_so_far] = []
        if self.temporal_phase_resolution > beat_time:
            self.beats[self.n_beats_so_far].append(self.time_so_far+beat_time)
        while time_within_beat < beat_time:
            if time_within_beat + self.temporal_phase_resolution > beat_time:
                self.beats[self.n_beats_so_far].append(self.time_so_far+beat_time)
                break
            else:
                time_within_beat += self.temporal_phase_resolution
                self.beats[self.n_beats_so_far].append(self.time_so_far+time_within_beat)
