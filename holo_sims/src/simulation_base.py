import holoocean
# import matplotlib.pyplot as plt
import numpy as np

import pyqtgraph as pqg
from pyqtgraph.Qt import QtCore


class ControllerBase:

    def __init__(self):
        self.name = "null"
        self.state = None

    def update_state(self, state):
        self.state = state

    @staticmethod
    def update_command():
        return np.array([0, 0, 0, 0, 50])


def simulate(scenario="OpenWater-Torpedo", controller=None, title="null"):
    # ==== CONFIG ====
    config = holoocean.packagemanager.get_scenario(scenario)
    # sensor_config = config['agents'][0]['sensors'][-1]["configuration"]

    # Set up
    print("Starting Simulation!")
    sensors_description = "".join([f"\n\t\t{sensor['sensor_type']}" for sensor in config['agents'][0]['sensors']])
    print(f"\tName: {config['name']}\n\tWorld: {config['world']}\n\tAgent: {config['agents'][0]['agent_type']}\n\tticks/s: {config['ticks_per_sec']}"
          f"\n\tSensors:{sensors_description}\n\tStarting location: {config['agents'][0]['location']}\n\tStarting Rotation: {config['agents'][0]['rotation']}")
    if controller is not None:
        print(f"Controller {controller.name} loaded!")
    else:
        print(f"No controller defined, using null controller")
        controller = ControllerBase()

    title = f"DSE07 | SEAMOS SIMULATION | {title} | {controller.name}"

    def rot2eul(R):
        beta = -np.arcsin(R[2,0])
        alpha = np.arctan2(R[2,1]/np.cos(beta),R[2,2]/np.cos(beta))
        gamma = np.arctan2(R[1,0]/np.cos(beta),R[0,0]/np.cos(beta))
        return np.array((alpha, beta, gamma))


    # pyQTgraph init
    app = pqg.mkQApp(title)
    plot_window = pqg.GraphicsLayoutWidget(show=True, title=title)
    plot_window.resize(1000, 600)
    plot_window.setWindowTitle(title)
    pqg.setConfigOptions(antialias=True)

    input_titles = ["right", "top", "left", "bottom", "thruster"]
    truth_titles = ["pos z", "vel", "roll", "pitch", "yaw"]
    imu_titles = ["accel x", "accel y", "accel z", "ang vel roll", "ang vel pitch", "ang vel yaw"]

    input_plots = []
    input_curves = []
    input_states = [[] for _ in range(len(input_titles))]

    truth_plots = []
    truth_curves = []
    truth_xy_states = [[], []]
    truth_states = [[] for _ in range(len(truth_titles))]

    imu_plots = []
    imu_curves = []
    imu_states = [[] for _ in range(len(imu_titles))]

    for plot_i in range(len(input_titles)):
        input_plots.append(plot_window.addPlot(title="input: " + input_titles[plot_i]))
        input_curves.append(input_plots[plot_i].plot())

    plot_window.nextRow()

    truth_xy_plot = plot_window.addPlot(title="truth: pos xy")
    truth_xy_curve = truth_xy_plot.plot()
    for plot_i in range(len(truth_titles)):
        truth_plots.append(plot_window.addPlot(title="truth: " + truth_titles[plot_i]))
        truth_curves.append(truth_plots[plot_i].plot())

    plot_window.nextRow()

    for plot_i in range(len(imu_titles)):
        imu_plots.append(plot_window.addPlot(title="imu: " + imu_titles[plot_i]))
        imu_curves.append(imu_plots[plot_i].plot())

    # RUN SIMULATION
    with holoocean.make(scenario) as env:
        def update_sim():
            command = controller.update_command()

            env.act("auv0", command)
            state = env.tick()

            controller.update_state(state)

            # print(state)

            for plot_i in range(len(command)):
                input_states[plot_i].append(command[plot_i])
                input_curves[plot_i].setData(input_states[plot_i])

            truth_state = state['PoseSensor']

            truth_euler = rot2eul(truth_state[0:3, 0:3])
            truth_xy_states[0].append(truth_state[0, 3])
            truth_xy_states[1].append(truth_state[1, 3])
            truth_states[0].append(truth_state[2, 3])

            velocities = state['VelocitySensor']
            truth_states[1].append((velocities[0]**2 + velocities[1]**2 + velocities[2]**2)**0.5)

            truth_states[2].append(truth_euler[0])
            truth_states[3].append(truth_euler[1])
            truth_states[4].append(truth_euler[2])

            truth_xy_curve.setData(x=truth_xy_states[0], y=truth_xy_states[1])
            for plot_i in range(len(truth_states)):
                truth_curves[plot_i].setData(truth_states[plot_i])

            imu_state = state['IMUSensor']
            imu_states[0].append(imu_state[0, 0])
            imu_states[1].append(imu_state[0, 1])
            imu_states[2].append(imu_state[0, 2])

            imu_states[3].append(imu_state[1, 0])
            imu_states[4].append(imu_state[1, 1])
            imu_states[5].append(imu_state[1, 2])

            for plot_i in range(len(imu_states)):
                imu_curves[plot_i].setData(imu_states[plot_i])

        timer = QtCore.QTimer()
        timer.timeout.connect(update_sim)
        timer.start()

        pqg.exec()

    print("Finished Simulation!")


if __name__ == "__main__":
    simulate()
