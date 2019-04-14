import numpy as np
import utility.utilmath as um
import utility.utilfile as uf
import utility.utilldg as ldg


class LunarDescentGuidance(object):

    def __init__(self, mission):
        self.inputs = uf.read_list('files/input_files/lunar_descent_guidance_inputs.txt')
        targets = uf.read_list('files/output_files/lunar_descent_guidance_targets.txt')

        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.planet = mission.conn.space_center.bodies[mission.planet]
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.ref_frame = self.moon.non_rotating_reference_frame

        # initial conditions
        self.latitude = self.inputs[9]
        self.longitude = self.inputs[10]
        self.landing_position_initial = (self.moon.equatorial_radius + self.inputs[11]) / self.moon.equatorial_radius * np.array(self.moon.surface_position(self.latitude, self.longitude, self.ref_frame))
        self.position_initial = np.array(self.vessel.position(self.ref_frame))
        self.velocity_initial = np.array(self.vessel.velocity(self.ref_frame))
        self.time_initial_ut = self.conn.space_center.ut

        self.guide_time_initial = -500
        self.time_initial = 0

        # state
        self.guide_time = None
        self.time = None


        # simulation parameters
        self.k = 0

        # self.vessel parameters
        self.leadtime = self.inputs[14]
        self.aftrim = self.inputs[16]

        # body parameters
        self.gravitational_parameter = self.inputs[4]
        self.surface_gravity = self.inputs[5]
        self.body_rotation = [0, -self.inputs[6], 0]

        self.kgt = [self.inputs[22], self.inputs[23], self.inputs[24], self.inputs[25]]

        # ullage burn time, trim burn time
        self.tullage = self.inputs[18]
        self.ttrim = self.inputs[17]
        
        self.tapi = self.inputs[26]
        self.tbrf = self.inputs[39]
        self.tapf = self.inputs[28]

        # approch phase states: initial position, midpoint position, midpoint velocity, terminal altitude, terminal altitude rate
        self.rapig = [self.inputs[29], 0, self.inputs[30]]
        self.rapmg = [self.inputs[31], 0, self.inputs[32]]
        self.vapmg = [self.inputs[33], 0, self.inputs[34]]
        self.rapfgx = self.inputs[35]
        self.vapfgx = self.inputs[36]

        # terminal constraints: terminal thrust angle, terminal thrust
        self.pbrf = self.inputs[41]
        self.fbrf = self.inputs[42]

        # guidance states
        self.brig, self.brtg, self.aptg = self.format_targets(targets)

    @staticmethod
    def format_targets(targets):
        brig = np.zeros(shape=(5, 3))
        brtg = np.zeros(shape=(5, 3))
        aptg = np.zeros(shape=(5, 3))

        for i in range(0, 5):
            for j in range(0, 3):
                brig[i][j] = targets[3 * i + j]
                brtg[i][j] = targets[3 * i + j + 15]
                aptg[i][j] = targets[3 * i + j + 30]

        return brig, brtg, aptg

    def execute(self):
        ldg.p63_ignition_algorithm(self)
        burn_time = self.time + self.time_initial_ut

        self.vessel.auto_pilot.engage()
        self.vessel.auto_pilot.reference_frame = self.ref_frame

        self.conn.space_center.rails_warp_factor = 4

        while self.conn.space_center.ut < burn_time - 90:
            self.ui.update_phase_current_time(int(abs(self.conn.space_center.ut - burn_time)))

        self.conn.space_center.rails_warp_factor = 2
        while self.conn.space_center.ut < burn_time - 30:
            self.ui.update_phase_current_time(int(abs(self.conn.space_center.ut - burn_time)))

        self.conn.space_center.rails_warp_factor = 0

        self.ui.update_phase_current('P63')
        self.ui.update_phase_next('P64')

        while True:
            position = np.array(self.vessel.position(self.ref_frame))
            velocity = np.array(self.vessel.velocity(self.ref_frame))
            landing_site = (self.moon.equatorial_radius + self.inputs[11]) / self.moon.equatorial_radius * np.array(self.moon.surface_position(self.latitude, self.longitude, self.ref_frame))
            afcp, tg, rg, vg = ldg.guide(self.guide_time, 0, position, velocity, landing_site, self)
            self.vessel.auto_pilot.target_direction = afcp
            maxafcp = self.vessel.max_thrust / self.vessel.mass
            ts = um.magnitude(afcp) / maxafcp
            if self.conn.space_center.ut > burn_time:
                if ts > 1:
                    ts = 1
                self.vessel.control.throttle = ts
                self.ui.update_phase_current_time(int(abs(tg - self.tbrf)))
            else:
                self.ui.update_phase_current_time(int(abs(self.conn.space_center.ut - burn_time)))
            if tg > self.tbrf:
                break

        print('swap to approach phase')
        self.ui.update_phase_current('P64')
        self.ui.update_phase_next('P66')
        self.target = self.aptg
        self.guidetime = self.tapi
        while tg < self.tapf and self.vessel.flight().surface_altitude > 30:
            self.ui.update_phase_current_time(int(abs(tg - self.tapf)))
            position = np.array(self.vessel.position(self.ref_frame))
            velocity = np.array(self.vessel.velocity(self.ref_frame))
            landing_site = (self.moon.equatorial_radius + self.inputs[11]) / self.moon.equatorial_radius * np.array(self.moon.surface_position(self.latitude, self.longitude, self.ref_frame))
            afcp, tg, rg, vg = ldg.guide(self.guide_time, 0, position, velocity, landing_site, self)
            self.vessel.auto_pilot.target_direction = afcp
            maxafcp = self.vessel.max_thrust / self.vessel.mass
            self.vessel.control.throttle = um.magnitude(afcp) / maxafcp

        self.ui.update_phase_current('P88')
        self.ui.update_phase_next('N/A')
        self.ui.update_phase_current_time(0)
        self.vessel.control.throttle = 0.00

        ref_frame = self.moon.reference_frame

        gm = self.moon.gravitational_parameter
        self.vessel.auto_pilot.engage()
        self.vessel.auto_pilot.reference_frame = ref_frame
        self.vessel.control.gear = True

        while self.vessel.flight().surface_altitude > 1.5:

            r = np.array(self.vessel.position(self.moon.reference_frame))
            vs = np.array(self.vessel.flight(ref_frame).vertical_speed)
            fa = np.array(self.vessel.direction(ref_frame))
            g = -gm / um.magnitude(r) ** 3 * r
            up = np.array(self.conn.space_center.transform_direction([1, 0, 0], self.vessel.surface_reference_frame, self.moon.reference_frame))
            v = np.array(self.vessel.velocity(self.moon.reference_frame))
            a = um.vector_exclude(up, -v)
            afyz = a / 10
            if self.vessel.flight().surface_altitude > 5:
                if um.magnitude(afyz) > 0.35 * um.magnitude(g):
                    afyz = 0.35 * g * um.normalize(afyz)
                if vs < 0:
                    afx = (((-1.00 - vs) / 1.5) * up - g) / np.cos(um.vector_angle(up, fa))
                else:
                    afx = 0
                afcp = afyz + afx
            else:
                if um.magnitude(afyz) > 0.35 * um.magnitude(g):
                    afyz = 0.35 * g * um.normalize(afyz)
                if vs < 0:
                    afx = (((-0.50 - vs) / 1.5) * up - g) / np.cos(um.vector_angle(up, fa))
                else:
                    afx = 0
                afcp = afyz + afx

            maxafcp = self.vessel.max_thrust / self.vessel.mass
            self.vessel.auto_pilot.target_direction = afcp
            self.vessel.control.throttle = um.magnitude(afcp) / maxafcp

        self.vessel.control.throttle = 0
        self.vessel.auto_pilot.disengage()



