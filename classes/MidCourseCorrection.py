import numpy as np
import subprocess
import utility.utilfile as uf
import utility.utilksp as uk


class MidCourseCorrection:

    def __init__(self, mission):
        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.planet = mission.conn.space_center.bodies[mission.planet]
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.ref_frame = self.planet.non_rotating_reference_frame

        # output parameters
        self.delta_v = None

    def calculate_correction(self, tli=False):
        moon_semi_major_axis = self.moon.orbit.semi_major_axis

        time_scale = (((self.planet.gravitational_parameter + self.moon.gravitational_parameter) / moon_semi_major_axis ** 3) ** -0.5)
        mcc_parameters = uf.read_list("files/output_files/mcc_parameters.txt")
        trajectory_solution = uf.read_list("files/output_files/trajectory_solution.txt")

        position_vessel = np.array(self.vessel.position(self.ref_frame)) / moon_semi_major_axis
        velocity_vessel = np.array(self.vessel.velocity(self.ref_frame)) / moon_semi_major_axis * time_scale

        position_moon = np.array(self.moon.position(self.ref_frame)) / moon_semi_major_axis
        velocity_moon = np.array(self.moon.velocity(self.ref_frame)) / moon_semi_major_axis * time_scale

        file = open("files/output_files/mid_course_correction_inputs.txt", "w")
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
        file.write(str(mcc_parameters[1] * self.planet.equatorial_radius / moon_semi_major_axis) + "\n")
        if tli:
            file.write(str((mcc_parameters[0] - self.conn.space_center.ut) / time_scale + trajectory_solution[6]) + "\n")
        else:
            file.write(str((mcc_parameters[0] - self.conn.space_center.ut) / time_scale) + "\n")
        file.close()

        p = subprocess.Popen([r"files/exe_files/mcc_km.exe"], cwd='files/output_files')
        p.wait()

        dv_data = np.array(uf.read_list("files/output_files/mid_course_correction_maneuver.txt")) * moon_semi_major_axis / time_scale
        dv_data = np.array([dv_data[0], dv_data[2], dv_data[1]])

        self.delta_v = dv_data - self.vessel.velocity(self.ref_frame)

    def execute_correction(self):
        uk.burn(self.conn, self.ui, self.ref_frame, self.delta_v, stopping_time=(1.0, 0.5, 1.0))

    def execute_correction_rcs(self):
        uk.burn(self.conn, self.ui, self.ref_frame, self.delta_v, stopping_time=(1.0, 0.5, 1.0), main_engine=False, specified_thrust=4000, specified_isp=240)




