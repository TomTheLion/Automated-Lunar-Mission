import numpy as np
import utility.utilmath as um
import utility.utilfile as uf
import utility.utilldg as ldg

class LunarDescentSimulation(object):

    def __init__(self, file):
        inputs = uf.read_list(file)

        # initial conditions
        self.landing_position_initial = ldg.landing_site(inputs[3] + inputs[11], np.radians(inputs[9]), np.radians(inputs[10]))

        self.guide_time_initial = -500
        self.time_initial = 0

        self.position_initial, self.velocity_initial = ldg.state_from_orbit(inputs[7], inputs[8], np.radians(inputs[9]), np.radians(inputs[10]), inputs[12], inputs[13], inputs[4])
        self.mass_initial = inputs[0]

        # state
        self.guide_time = None
        self.time = None

        self.position = None
        self.velocity = None
        self.mass = None

        # simulation parameters
        self.k = 0
        self.throttflag = False
        self.thrott = None
        self.throtta = None
        self.mbrf = inputs[40]
        self.jbrfg = [0, 0, 0]
        self.sbrfg = [0, 0, 0]

        # vessel parameters
        self.leadtime = inputs[14]
        self.maxthrust = inputs[1]
        self.vexh = inputs[2]
        self.minthrottle = inputs[15]
        self.aftrim = inputs[16]

        # body parameters
        self.gravitational_parameter = inputs[4]
        self.surface_gravity = inputs[5]
        self.body_rotation = [0, 0, inputs[6]]

        # guidance constants: transport delay, jerk coefficient, ignition coefficient, guide time coefficient
        self.tau = inputs[19]
        self.kj = inputs[20]
        self.kig = inputs[21]
        self.kgt = [inputs[22], inputs[23], inputs[24], inputs[25]]

        # ullage burn time, trim burn time
        self.tullage = inputs[18]
        self.ttrim = inputs[17]

        # descent guidance parameters

        # braking phase guide times: initial, throttle recovery, final
        self.tbri = inputs[37]
        self.thrott = inputs[38]
        self.tbrf = inputs[39]

        # approach phase guide times: initial, midpoint, final
        self.tapi = inputs[26]
        self.tapm = inputs[27]
        self.tapf = inputs[28]

        # approch phase states: initial position, midpoint position, midpoint velocity, terminal altitude, terminal altitude rate
        self.rapig = [inputs[29], 0, inputs[30]]
        self.rapmg = [inputs[31], 0, inputs[32]]
        self.vapmg = [inputs[33], 0, inputs[34]]
        self.rapfgx = inputs[35]
        self.vapfgx = inputs[36]

        # terminal constraints: terminal thrust angle, terminal thrust
        self.pbrf = inputs[41]
        self.fbrf = inputs[42]

        # guidance states
        self.apig = None  # approach phase ignition state
        self.aptg = None  # approach phase target state
        self.brig = None  # braking phase initial state
        self.brtg = None  # braking phase target state
        self.brfg = None  # braking phase final state


    # simulation
    def solve(self):
        while True:
            # set apig, aptg
            ldg.approach_phase_targeting_routine(self)
            # set brig, brtg
            ldg.breaking_phase_targeting_routine(self)
            # # set time, guide_time
            ldg.p63_ignition_algorithm(self)

            # set position and velocity at ignition
            self.time += self.tullage
            self.position, self.velocity = um.cse_rv(self.position_initial, self.velocity_initial, self.time, self.gravitational_parameter)
            self.mass = self.mass_initial
            # # step forward, delta time is trim time

            ldg.simulate(self, self.brtg, self.ttrim, 1, self.aftrim, False)
            ldg.simulate(self, self.brtg, self.tbrf, 1, -1, True)
            conv = self.brf_update()
            self.reset()
            if conv:
                self.set_brig()
                uf.write_list(np.array([self.brig, self.brtg, self.aptg]).flatten(), 'files/output_files/lunar_descent_guidance_targets.txt')
                break

    def brf_update(self):
        lt = um.vector_axis_angle(self.landing_position_initial, np.array([0, 0, 1]), self.body_rotation[2] * self.time)
        afcp, self.guide_time, rg, vg = ldg.guide(self.guide_time, 0, self.position, self.velocity, lt, self)
        self.brig[0][2] += self.kig * (self.throtta - self.thrott)
        self.mbrf = self.mass * np.exp((um.magnitude(self.brfg[1]) - um.magnitude(vg)) / self.vexh)
        jbrtga, sbrtga = ldg.eq2109(self.guide_time, [self.brtg[0], self.brtg[1], self.brtg[2], rg, vg])
        brfga = np.matmul(um.state_transition_matrix(self.tbrf, 0), [self.brtg[0], self.brtg[1], self.brtg[2], [jbrtga[0], jbrtga[1], self.brtg[3][2]], sbrtga])
        self.jbrfg = brfga[3]
        self.sbrfg = brfga[4]

        q1 = abs(brfga[2][2] - self.brfg[2][2]) / um.magnitude(self.brfg[2])
        q2 = abs(brfga[2][0] - self.brfg[2][0]) / um.magnitude(self.brfg[2])
        q3 = abs(brfga[1][0] - self.brfg[1][0]) / um.magnitude(self.brfg[1])
        q4 = abs(self.throtta - self.thrott)

        print(q1, q2, q3, q4)

        if q1 > 1e-4:
            return False
        if q2 > 1e-4:
            return False
        if q3 > 1e-4:
            return False
        if q4 > 1e-1:
            return False

        return True

    def reset(self):
        self.guide_time = self.guide_time_initial
        self.time = self.time_initial

        self.position = self.position_initial
        self.velocity = self.velocity_initial
        self.mass = self.mass_initial

        self.throttflag = False
        return

    def set_brig(self):
        ldg.p63_ignition_algorithm(self)
        self.time += self.tullage + self.ttrim
        self.position, self.velocity = um.cse_rv(self.position_initial, self.velocity_initial, self.time, self.gravitational_parameter)
        self.target = self.brtg

        lt = um.vector_axis_angle(self.landing_position_initial, np.array([0, 0, 1]), self.body_rotation[2] * self.time)
        afcp, self.guide_time, rg, vg = ldg.guide(self.guide_time, 0, self.position, self.velocity, lt, self)
        self.brig[0][0] = rg[0]
        self.brig[1] = vg

        self.reset()
