"""
Simple PID implementation (by Lilly)

Not verified
"""

class PID:

    def __init__(self, p_gain, i_gain, d_gain, value_range):
        self.p_gain = p_gain
        self.i_gain = i_gain
        self.d_gain = d_gain
        self.set_point = 0
        self.value_range = value_range

        self.previous_error = 0
        self.i_value = 0

    def set(self, set_point):
        self.set_point = set_point

    def get_output(self, current_value):
        error = self.set_point - current_value
        self.i_value += error
        self.i_value = min(max(self.i_value, self.value_range[0]), self.value_range[1])
        d_value = error - self.previous_error
        self.previous_error = error
        return min(max(error * self.p_gain + self.i_value * self.i_gain + d_value * self.d_gain, self.value_range[0]), self.value_range[1])
