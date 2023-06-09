"""
Example implementation of sim based on simulation_base

Has a PID controller for the yaw and a state estimator for roll, pitch and yaw
"""


from simulation_base import simulate, ControllerBase, rot2eul, wrap_around
import pyqtgraph as pqg
import numpy as np
from simple_pid import PID


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
