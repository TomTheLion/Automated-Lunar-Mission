import numpy as np
import utility.utilfile as uf
import utility.utilmath as um
import utility.utilksp as uk


class LunarOrbitInsertion:

    def __init__(self, mission):
        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.ref_frame = self.moon.non_rotating_reference_frame

        self.normal = np.cross(self.vessel.position(self.ref_frame), self.vessel.velocity(self.ref_frame))
        self.normal /= np.linalg.norm(self.normal)
        self.inclination = self.vessel.orbit.inclination

        self.latitude, self.longitude, self.apoapsis, self.periapsis = uf.read_list('files/input_files/lunar_landing_inputs.txt')
        self.lan = self.longitude

    def circularize(self):
        self.conn.space_center.warp_to(self.conn.space_center.ut + self.vessel.orbit.time_to_periapsis - 15)

        circular_speed = np.sqrt(self.moon.gravitational_parameter / np.linalg.norm(np.array(self.vessel.position(self.ref_frame))))
        circular_direction = um.normalize(np.cross(np.cross(np.array(self.vessel.position(self.ref_frame)), np.array(self.vessel.velocity(self.ref_frame))), np.array(self.vessel.position(self.ref_frame))))
        delta_velocity = circular_speed * circular_direction - np.array(self.vessel.velocity(self.ref_frame))

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(1.0, 0.5, 1.0))

    def change_inclination(self):
        sma_1 = (self.vessel.orbit.semi_major_axis + self.apoapsis + self.moon.equatorial_radius) / 2
        sma_2 = self.apoapsis + self.moon.equatorial_radius
        sma_3 = (self.apoapsis + self.periapsis) / 2 + self.moon.equatorial_radius
        period_1 = 2 * np.pi * (sma_1 ** 3 / self.moon.gravitational_parameter) ** 0.5
        period_2 = 2 * np.pi * (sma_2 ** 3 / self.moon.gravitational_parameter) ** 0.5
        period_3 = 2 * np.pi * (sma_3 ** 3 / self.moon.gravitational_parameter) ** 0.5

        self.lan = self.longitude + (period_1 / 2 + period_2 + period_3 / 2) / self.moon.orbit.period * 360 + np.degrees(self.moon.rotation_angle)

        delta_longitude = um.safe_asin(np.tan(np.radians(self.latitude)) / np.tan(self.inclination))

        closest_normal = self.calculate_landing_normal(np.radians(self.lan), delta_longitude)

        angle_to_maneuver = self.calculate_angle_to_maneuver(closest_normal)

        self.conn.space_center.warp_to(self.conn.space_center.ut + angle_to_maneuver / (2 * np.pi) * self.vessel.orbit.period - 15)

        attitude = um.normalize(np.cross(closest_normal, np.array(self.vessel.position(self.ref_frame))))
        delta_velocity = attitude * np.sqrt(self.moon.gravitational_parameter * (2 / np.linalg.norm(self.vessel.position(self.ref_frame)) - 1 / sma_1)) - self.vessel.velocity(self.ref_frame)

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(1.0, 0.5, 1.0))

    def calculate_landing_normal(self, longitude_of_ascending_node, delta_longitude):
        ref_axis = -np.array(self.vessel.orbit.reference_plane_normal(self.ref_frame))
        ref_direction = np.array(self.vessel.orbit.reference_plane_direction(self.ref_frame))
        node_1 = um.vector_axis_angle(ref_direction, ref_axis, longitude_of_ascending_node + delta_longitude)
        node_2 = um.vector_axis_angle(ref_direction, ref_axis, longitude_of_ascending_node - delta_longitude)
        velocity_1 = um.vector_axis_angle(np.cross(ref_axis, node_1), node_1, self.inclination)
        velocity_2 = um.vector_axis_angle(np.cross(ref_axis, node_2), node_2, -self.inclination)

        normal_1p = np.cross(node_1, velocity_1)
        normal_1n = np.cross(node_1, -velocity_1)
        normal_2p = np.cross(node_2, velocity_2)
        normal_2n = np.cross(node_2, -velocity_2)

        angle_1p = um.vector_angle(normal_1p, self.normal)
        angle_1n = um.vector_angle(normal_1n, self.normal)
        angle_2p = um.vector_angle(normal_2p, self.normal)
        angle_2n = um.vector_angle(normal_2n, self.normal)

        angles = [angle_1p, angle_1n, angle_2p, angle_2n]
        normals = [normal_1p, normal_1n, normal_2p, normal_2n]

        return normals[angles.index(min(angles))]

    def calculate_angle_to_maneuver(self, normal):
        node_1 = np.cross(self.normal, normal)
        node_2 = -node_1

        vector_a_1 = np.cross(node_1, self.vessel.position(self.ref_frame))
        vector_a_2 = np.cross(node_2, self.vessel.position(self.ref_frame))
        vector_b = np.cross(self.vessel.position(self.ref_frame), self.vessel.velocity(self.ref_frame))

        angle_to_1 = um.vector_angle(node_1, self.vessel.position(self.ref_frame))
        angle_to_2 = um.vector_angle(node_2, self.vessel.position(self.ref_frame))

        if np.dot(vector_a_1, vector_b) > 0:
            angle_to_1 = 2 * np.pi - angle_to_1

        if np.dot(vector_a_2, vector_b) > 0:
            angle_to_2 = 2 * np.pi - angle_to_2

        if angle_to_1 < angle_to_2:
            print('time_1 is closer')
            angle_to_maneuver = angle_to_1
        else:
            print('time_2 is closer')
            angle_to_maneuver = angle_to_2

        return angle_to_maneuver

    def calculate_descent_orbit_time(self):
        ref_axis = np.array(self.vessel.orbit.reference_plane_normal(self.ref_frame))
        normal = um.normalize(np.cross(self.vessel.position(self.ref_frame), self.vessel.velocity(self.ref_frame)))
        node = np.cross(normal, ref_axis)
        landing_site = self.moon.surface_position(self.latitude, self.longitude, self.ref_frame)

        angle_dir = 1
        if np.dot(node, landing_site) > 0:
            angle_dir = -1
            node = -node

        theta = um.safe_asin(np.sin(np.radians(self.latitude)) / np.sin(self.inclination))
        target = um.vector_axis_angle(node, normal, theta * angle_dir)
        target_angle = um.vector_angle(self.vessel.position(self.ref_frame), target)
        delta_direction = np.dot(normal, np.cross(target, self.vessel.position(self.ref_frame)))

        if delta_direction > 0:
            target_angle = 2 * np.pi - target_angle

        target_time = target_angle * self.vessel.orbit.period / (2 * np.pi)
        return target_time

    def descent_orbit(self):
        target_time = self.calculate_descent_orbit_time()
        self.conn.space_center.warp_to(self.conn.space_center.ut + target_time - 15)

        position = np.array(self.vessel.position(self.ref_frame))
        velocity = np.array(self.vessel.velocity(self.ref_frame))

        speed = np.sqrt(self.moon.gravitational_parameter * (2 / np.linalg.norm(position) - 2 / (np.linalg.norm(position) + self.moon.equatorial_radius + self.periapsis)))
        direction = um.normalize(np.cross(np.cross(np.array(self.vessel.position(self.ref_frame)), np.array(self.vessel.velocity(self.ref_frame))), np.array(self.vessel.position(self.ref_frame))))

        delta_velocity = speed * direction - velocity

        uk.burn(self.conn, self.ui, self.ref_frame, delta_velocity, stopping_time=(0.1, 0.5, 0.1))
