from simulation_base import simulate, ControllerBase
import pyqtgraph as pqg
import numpy as np


def rot2eul(R):
    beta = -np.arcsin(R[2, 0])
    alpha = np.arctan2(R[2, 1] / np.cos(beta), R[2, 2] / np.cos(beta))
    gamma = np.arctan2(R[1, 0] / np.cos(beta), R[0, 0] / np.cos(beta))
    return np.array((alpha, beta, gamma)) / np.pi * 180


def wrap_around(value, bounds):
    difference = bounds[1] - bounds[0]
    if bounds[0] <= value <= bounds[1]:
        return value
    elif value < bounds[0]:
        return value + difference
    elif value > bounds[1]:
        return value - difference
    else:
        return value


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


class IdiotIMUController(ControllerBase):

    def __init__(self):
        super().__init__()
        self.name = "idiot_imu_estimator"
        self.imu_states = []
        self.sim_i = 0

        self.state_estimation = [0, 0, 0, 0, 0, 0]
        self.states = [[] for _ in range(len(self.state_estimation))]
        self.err_states = [[] for _ in range(len(self.state_estimation))]

        self.has_plots = True
        self.plots = []
        self.curves = []
        self.err_plots = []
        self.err_curves = []

        self.roll_pid = PID(0.1, 0.05, 0, (-20, 20))
        self.roll_pid.set(0)

        self.yaw_pid = PID(0.1, 0.05, -0.1, (-20, 20))
        self.yaw_pid.set(90)

        self.yaw_mix = np.array([0, 0, 0, 0, 1, -1, -1, 1])

        self.state_truth = [0, 0, 0, 0, 0, 0]

    def init_plots(self, plot_window):
        titles = ["pos x", "pos y", "pos z", "roll", "pitch", "yaw"]
        for plot_i in range(len(self.states)):
            self.plots.append(plot_window.addPlot(title="controller: SE " + titles[plot_i]))
            self.curves.append(self.plots[plot_i].plot())
        plot_window.nextRow()
        for plot_i in range(len(self.states)):
            self.err_plots.append(plot_window.addPlot(title="controller: SE err " + titles[plot_i]))
            self.err_curves.append(self.err_plots[plot_i].plot(pen='r'))

    def update_state(self, state):
        self.sim_i += 1
        # self.imu_states.append(state['IMUSensor'])
        truth_state = state['PoseSensor']
        self.state_truth = np.concatenate((truth_state[0:3, 3].T, rot2eul(truth_state[0:3, 0:3])))
        imu_state = state['IMUSensor']
        sensor_frametime = 1/200

        for plot_i in [3, 4, 5]:
            self.state_estimation[plot_i] += sensor_frametime * - imu_state[1, plot_i - 3] / np.pi * 180
            self.state_estimation[plot_i] = wrap_around(self.state_estimation[plot_i], (-180, 180))

        for plot_i in range(len(self.state_estimation)):
            self.states[plot_i].append(self.state_estimation[plot_i])
            self.curves[plot_i].setData(self.states[plot_i])
            self.err_states[plot_i].append(self.state_estimation[plot_i] - self.state_truth[plot_i])
            self.err_curves[plot_i].setData(self.err_states[plot_i])

    def update_command(self):
        yaw_command = self.yaw_pid.get_output(self.state_truth[5])
        return yaw_command * self.yaw_mix + np.array([0, 0, 0, 0, 100, 100, 100, 100])

# 50 * np.sin(self.sim_i / 100) + 50]


simulate(controller=IdiotIMUController(), title='Idiot State Estimation')
