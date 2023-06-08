from simulation_base import simulate, ControllerBase


class IdiotIMUController(ControllerBase):

    def __init__(self):
        super().__init__()
        self.name = "idiot_imu_estimator"

    def update_state(self, state):
        pass


simulate(controller=IdiotIMUController())
