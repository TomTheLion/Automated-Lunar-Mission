import numpy as np
import time
import subprocess
import utility.utilmath as um
import utility.utilksp as uk

class Rendezvous:

    def __init__(self, mission):
        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.planet = mission.conn.space_center.bodies[mission.planet]
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.target = self.conn.space_center.target_vessel
        self.ref_frame = self.moon.non_rotating_reference_frame
        self.target_vessel = self.conn.space_center.target_vessel

    def find_orbit(self, p, plane_change_flag):

        r = self.vessel.position(self.ref_frame)
        v = self.vessel.velocity(self.ref_frame)

        rt = self.conn.space_center.target_vessel.position(self.ref_frame)
        vt = self.conn.space_center.target_vessel.velocity(self.ref_frame)

        rt2, vt2 = um.cse_rv(rt, vt, p, self.moon.gravitational_parameter)

        if plane_change_flag:
            v1 = np.sqrt(self.moon.gravitational_parameter * (2 / um.magnitude(r) - 2 / (um.magnitude(r) + um.magnitude(rt))))

            v1 *= um.normalize(np.cross(np.cross(r, v), r))
        else:
            v1, v2 = um.cse_rr(r, rt2, p, -1, self.moon.gravitational_parameter)

        delta_velocity = v1 - v

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(0.1, 0.5, 0.1))

    def match_speed(self, d, v):
        while um.magnitude(self.target.position(self.vessel.reference_frame)) > d:
            pass

        cv = np.array(self.vessel.velocity(self.ref_frame)) - np.array(self.target.velocity(self.ref_frame))
        cp = np.array(self.vessel.position(self.ref_frame)) - np.array(self.target.position(self.ref_frame))
        desv = um.normalize(-cp) * v

        delta_velocity = desv - cv

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(0.25, 0.5, 0.25))

    def match_speed_mono(self, d, v):
        while um.magnitude(self.target.position(self.vessel.reference_frame)) > d:
            pass

        cv = np.array(self.vessel.velocity(self.ref_frame)) - np.array(self.target.velocity(self.ref_frame))
        cp = np.array(self.vessel.position(self.ref_frame)) - np.array(self.target.position(self.ref_frame))
        desv = um.normalize(-cp) * v

        delta_velocity = desv - cv

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(0.25, 0.5, 0.25), main_engine=False, specified_thrust=4000, specified_isp=240, facing=-1)
        self.vessel_attitude()

    def execute_rendezvous(self):
        r = self.vessel.position(self.ref_frame)
        v = self.vessel.velocity(self.ref_frame)

        rt = self.target_vessel.position(self.ref_frame)

        ra = self.target_vessel.orbit.semi_major_axis
        rp = self.vessel.orbit.semi_major_axis

        p1 = self.target_vessel.orbit.period
        p2 = 2 * np.pi * np.sqrt(((ra + rp) / 2) ** 3 / self.moon.gravitational_parameter)
        p3 = self.vessel.orbit.period

        transfer_time = p2 / 2
        target_transfer_angle = transfer_time / p1 * np.pi * 2

        angle_rate = 2 * np.pi / p3 - 2 * np.pi / p1

        angle_current = um.vector_angle(r, rt)

        angle_direction = np.dot(np.cross(r, v), np.cross(r, rt))

        if angle_current < np.pi - target_transfer_angle:
            if angle_direction > 0:
                delta_angle = angle_current + target_transfer_angle + np.pi
            else:
                delta_angle = -angle_current + target_transfer_angle + np.pi
        else:
            if angle_direction > 0:
                delta_angle = angle_current - (np.pi - target_transfer_angle)
            else:
                delta_angle = -angle_current + target_transfer_angle + np.pi

        tt = delta_angle / angle_rate

        self.conn.space_center.warp_to(self.conn.space_center.ut + tt)

        self.find_orbit(p2 / 2, True)

        self.conn.space_center.warp_to(self.conn.space_center.ut + p2 / 2 - 120)

        self.match_speed(10000, 50)
        self.match_speed(2000, 30)
        self.match_speed(1000, 15)

        self.vessel.control.rcs = True

        self.vessel.auto_pilot.engage()
        self.vessel.auto_pilot.reference_frame = self.ref_frame
        self.target_vessel.auto_pilot.engage()
        self.target_vessel.auto_pilot.reference_frame = self.ref_frame
        self.vessel.auto_pilot.stopping_time = (0.25, 0.5, 0.25)
        self.vessel_attitude()

        self.match_speed_mono(500, 10)
        self.match_speed_mono(200, 3)
        self.match_speed_mono(65, 1)
        self.match_speed_mono(25, 0.5)
        self.match_speed_mono(12, 0.2)

    def vessel_attitude(self):
        target_direction = np.array(self.target.position(self.ref_frame)) - np.array(self.vessel.position(self.ref_frame))
        self.vessel.auto_pilot.engage()
        self.vessel.auto_pilot.reference_frame = self.ref_frame
        self.target_vessel.auto_pilot.engage()
        self.target_vessel.auto_pilot.reference_frame = self.ref_frame
        self.vessel.auto_pilot.target_direction = target_direction
        self.target_vessel.auto_pilot.target_direction = -target_direction
