import numpy as np
import time
import subprocess
import utility.utilfile as uf
import utility.utilksp as uk


class LunarReturn:

    def __init__(self, mission):
        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.planet = mission.conn.space_center.bodies[mission.planet]
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.ref_frame = self.planet.non_rotating_reference_frame

        # output parameters
        self.maneuver_velocity = None
        self.maneuver_time = None

    def calculate_return(self):
        moon_semi_major_axis = self.moon.orbit.semi_major_axis

        time_scale = (((self.planet.gravitational_parameter + self.moon.gravitational_parameter) / moon_semi_major_axis ** 3) ** -0.5)
        return_parameters = uf.read_list("files/input_files/return_parameters.txt")

        position_vessel = np.array(self.vessel.position(self.ref_frame)) / moon_semi_major_axis
        velocity_vessel = np.array(self.vessel.velocity(self.ref_frame)) / moon_semi_major_axis * time_scale

        position_moon = np.array(self.moon.position(self.ref_frame)) / moon_semi_major_axis
        velocity_moon = np.array(self.moon.velocity(self.ref_frame)) / moon_semi_major_axis * time_scale

        start_time = self.conn.space_center.ut

        file = open("files/output_files/lunar_return_inputs.txt", "w")
        file.write(str(position_vessel[0]) + "\n")
        file.write(str(position_vessel[2]) + "\n")
        file.write(str(position_vessel[1]) + "\n")
        file.write(str(velocity_vessel[0]) + "\n")
        file.write(str(velocity_vessel[2]) + "\n")
        file.write(str(velocity_vessel[1]) + "\n")
        file.write(str(position_moon[0]) + "\n")
        file.write(str(position_moon[2]) + "\n")
        file.write(str(position_moon[1]) + "\n")
        file.write(str(velocity_moon[0]) + "\n")
        file.write(str(velocity_moon[2]) + "\n")
        file.write(str(velocity_moon[1]) + "\n")
        file.write(str(2.0 * return_parameters[0] * self.planet.equatorial_radius / moon_semi_major_axis) + "\n")
        file.write(str(self.vessel.orbit.period / time_scale) + "\n")
        file.close()

        p = subprocess.Popen([r"files/exe_files/return_km.exe"], cwd='files/output_files')
        p.wait()

        return_data = np.array(uf.read_list("files/output_files/lunar_return_maneuver.txt"))
        self.maneuver_velocity = np.array([return_data[0], return_data[2], return_data[1]]) * moon_semi_major_axis / time_scale
        self.maneuver_time = return_data[3] * time_scale + start_time
        uf.write_list([start_time + return_data[4] * time_scale, return_parameters[0]], 'files/output_files/mcc_parameters.txt')

    def execute_correction(self):
        self.conn.space_center.warp_to(self.maneuver_time - 15)
        delta_v = self.maneuver_velocity - np.array(self.vessel.velocity(self.ref_frame))
        uk.burn(self.conn, self.ui, self.ref_frame, delta_v, stopping_time=(1.0, 0.5, 1.0))
