import numpy as np
import subprocess
import time
import utility.utilmath as um
import utility.utilfile as uf
import utility.utilksp as uk


class TransLunarInjection:

    def __init__(self, mission):
        # read inputs
        inputs = uf.read_list('files/input_files/trans_lunar_injection_inputs.txt')

        # krpc parameters
        self.conn = mission.conn
        self.ui = mission.ui
        self.vessel = mission.conn.space_center.active_vessel
        self.planet = mission.conn.space_center.bodies[mission.planet]
        self.moon = mission.conn.space_center.bodies[mission.moon]
        self.ref_frame = self.planet.non_rotating_reference_frame

        # control parameters
        self.fail = False
        self.convergence = False

        # vessel properties
        self.thrust = inputs[0]
        self.specific_impulse = inputs[1]
        self.initial_mass = inputs[2]
        self.launch_time_offset = inputs[3]
        self.launch_angle_offset = inputs[4]

        # trajectory constraints
        self.radius_initial = inputs[5]
        self.radius_leave = self.radius_initial
        self.radius_return = inputs[6]
        self.radius_approach = inputs[7]
        self.inclination_leave = inputs[8]
        self.inclination_return = inputs[9]
        self.launch_direction = inputs[10]

        # estimate time to approach
        self.time_approach = 0.5 * self.planet.rotational_period + 0.2 * self.moon.orbit.period

        # output variables
        self.vessel_orbital_parameters = None
        self.moon_orbital_parameters = None

        # launch and burn parameters
        self.time_launch = 0
        self.time_reach_orbit = None
        self.time_start_burn = None
        self.time_end_burn = None
        self.burn_end_radius = None
        self.burn_end_velocity = None
        self.burn_end_acceleration = None
        self.angle_burn = None
        self.time_start = None
        self.burn_delta_velocity = None
        self.update_constant = self.planet.rotational_period / (self.moon.orbit.period - self.planet.rotational_period)

    def solve_trajectory(self):
        while not self.convergence:
            self.vessel_orbital_parameters, self.moon_orbital_parameters, subprocess_return_code = self.calculate_trajectory()
            if subprocess_return_code != 0:
                self.fail = True
                break
            self.angle_burn, self.convergence = self.update_constraints()
            self.time_start = self.conn.space_center.ut
            uf.write_list([self.time_start + self.time_approach, self.radius_return], 'files/output_files/mcc_parameters.txt')

        if self.convergence:
            lan = np.degrees(self.moon.orbit.longitude_of_ascending_node) + self.vessel_orbital_parameters[3] - self.moon_orbital_parameters[3]
            if lan < 0:
                lan = 360 + lan

            self.ui.update_message('R: {:3.1f} I: {:3.1f} L: {:3.2f} D: {}'.format((self.radius_initial - 1) * self.planet.equatorial_radius / 1000, self.inclination_leave, lan, self.launch_direction[0]))

        print('Convergence: ', self.convergence)

    def calculate_trajectory(self):
        trajectory_constraints = []
        trajectory_constraints.append(self.radius_leave * self.planet.equatorial_radius / self.moon.orbit.semi_major_axis)
        trajectory_constraints.append(self.radius_return * self.planet.equatorial_radius / self.moon.orbit.semi_major_axis)
        trajectory_constraints.append(self.radius_approach * self.moon.equatorial_radius / self.moon.orbit.semi_major_axis)
        trajectory_constraints.append(np.radians(self.inclination_leave))
        trajectory_constraints.append(np.radians(self.inclination_return))
        trajectory_constraints.append(self.time_approach / (((self.planet.gravitational_parameter + self.moon.gravitational_parameter) / self.moon.orbit.semi_major_axis ** 3) ** -0.5))

        uf.write_list(trajectory_constraints, 'files/output_files/trajectory_constraints.txt')
        self.write_orbital_parameters(self.moon, 'files/output_files/lunar_orbital_parameters.txt')

        p = subprocess.Popen([r'files/exe_files/ur3bp_km.exe'], cwd='files/output_files')
        p.wait()

        if p.returncode == 0:
            return uf.read_list('files/output_files/orbital_parameters.txt'), uf.read_list('files/output_files/orbital_parameters_moon.txt'), p.returncode
        else:
            return 0, 0, p.returncode

    def update_constraints(self):
        initial_position = np.array([self.radius_initial * self.planet.equatorial_radius, 0, 0])
        initial_velocity = np.array([0, np.sqrt(self.planet.gravitational_parameter / initial_position[0]), 0])
        burn = self.Burn(initial_position, initial_velocity, self.initial_mass, self.thrust, self.specific_impulse, self.planet.gravitational_parameter)
        burn.integrate_to_eccentricity(self.vessel_orbital_parameters[1], 1)

        burn_sweep_angle = np.arctan2(burn.position[1], burn.position[0])
        burn_periapsis_sweep_angle = burn_sweep_angle - burn.true_anomaly()
        burn_angle_offset = self.vessel_orbital_parameters[4] * np.pi / 180 - burn_periapsis_sweep_angle
        delta_time_orbit = (burn_angle_offset - self.launch_angle_offset) * np.sqrt((self.radius_initial * self.planet.equatorial_radius) ** 3 / self.planet.gravitational_parameter)

        if delta_time_orbit < 0:
            delta_time_orbit += 2 * np.pi * np.sqrt((self.radius_initial * self.planet.equatorial_radius) ** 3 / self.planet.gravitational_parameter)

        inclination = self.vessel_orbital_parameters[2] * np.pi / 180

        longitude_of_ascending_node = (self.moon.orbit.longitude_of_ascending_node * 180 / np.pi + self.vessel_orbital_parameters[3] - self.moon_orbital_parameters[3]) * np.pi / 180

        delta_time_tli = (((self.planet.gravitational_parameter + self.moon.gravitational_parameter) / self.moon.orbit.semi_major_axis ** 3) ** -0.5) * self.vessel_orbital_parameters[5]
        delta_time_tli -= burn.mean_anomaly() * np.sqrt(burn.semi_major_axis() ** 3 / self.planet.gravitational_parameter)

        node_delta_north, node_delta_south = self.node_delta(inclination, longitude_of_ascending_node, self.vessel, self.ref_frame)
        if self.launch_direction == 'north':
            delta_time_launch = node_delta_north * self.planet.rotational_period / (2 * np.pi)
        elif self.launch_direction == 'south':
            delta_time_launch = node_delta_south * self.planet.rotational_period / (2 * np.pi)
        else:
            delta_time_launch = 0

        delta_time_burn = delta_time_launch + self.launch_time_offset + delta_time_orbit
        delta_time_approach = delta_time_burn + delta_time_tli

        if abs(delta_time_approach - self.time_approach) < 1 and abs(burn.periapsis() / self.planet.equatorial_radius - self.radius_leave) < 1e-5:
            convergence = True
        else:
            convergence = False

        if delta_time_launch < (self.time_approach - delta_time_approach) / self.update_constant:
            delta_time_approach += (1 + self.update_constant) * self.planet.rotational_period

        self.radius_leave = burn.periapsis() / self.planet.equatorial_radius
        self.time_launch = delta_time_launch
        self.time_reach_orbit = delta_time_launch + self.launch_time_offset
        self.time_start_burn = delta_time_burn
        self.time_end_burn = delta_time_burn + burn.time
        self.time_approach = delta_time_approach - self.update_constant * (self.time_approach - delta_time_approach)
        self.burn_end_radius = np.linalg.norm(burn.position)
        self.burn_end_velocity = np.linalg.norm(burn.velocity)
        self.burn_end_acceleration = burn.max_thrust / burn.mass
        self.burn_delta_velocity = np.log(burn.initial_mass / burn.mass) * burn.exhaust_velocity

        return burn_angle_offset, convergence

    @staticmethod
    def node_vector(inclination, vessel, ref_frame):
        position = np.array(vessel.position(ref_frame))
        latitude = np.radians(vessel.flight().latitude)
        ref_axis = np.array(vessel.orbit.reference_plane_normal(ref_frame))
        node_vector_projection_delta = um.safe_asin(np.tan(latitude) / np.tan(inclination))
        longitude_vector = um.vector_exclude(ref_axis, position)
        return um.vector_axis_angle(longitude_vector, ref_axis, node_vector_projection_delta), um.vector_axis_angle(longitude_vector, ref_axis, np.pi - node_vector_projection_delta)

    def node_delta(self, inclination, longitude_of_ascending_node, vessel, ref_frame):
        ref_axis = np.array(vessel.orbit.reference_plane_normal(ref_frame))
        ref_direction = np.array(vessel.orbit.reference_plane_direction(ref_frame))

        node_north, node_south = self.node_vector(inclination, vessel, ref_frame)
        target_node = um.vector_axis_angle(ref_direction, ref_axis, -longitude_of_ascending_node)

        node_delta_north = um.vector_angle(node_north, target_node)
        node_delta_south = um.vector_angle(node_south, target_node)
        delta_direction_north = np.dot(ref_axis, np.cross(target_node, node_north))
        delta_direction_south = np.dot(ref_axis, np.cross(target_node, node_south))

        if delta_direction_north < 0:
            node_delta_north = 2 * np.pi - node_delta_north
        if delta_direction_south < 0:
            node_delta_south = 2 * np.pi - node_delta_south

        return node_delta_north, node_delta_south

    @staticmethod
    def write_orbital_parameters(body, filename):
        file = open(filename, 'w')
        file.write(str(body.orbit.eccentricity) + '\n')
        file.write(str(body.orbit.inclination + 1e-12) + '\n')
        file.write(str(body.orbit.argument_of_periapsis) + '\n')
        file.write(str(body.orbit.true_anomaly) + '\n')
        file.close()

    class Burn:

        def __init__(self, position, velocity, mass, max_thrust, specific_impulse, gravitational_parameter):
            self.time = 0.0
            self.initial_position = position
            self.initial_velocity = velocity
            self.initial_mass = mass
            self.position = position
            self.velocity = velocity
            self.mass = mass
            self.attitude = np.cross(np.cross(position, velocity), position)
            self.attitude /= np.linalg.norm(self.attitude)
            self.max_thrust = max_thrust
            self.exhaust_velocity = specific_impulse * 9.80665
            self.gravitational_parameter = gravitational_parameter
            self.burn = True

        def reset(self):
            self.time = 0.0
            self.position = self.initial_position
            self.velocity = self.initial_velocity
            self.mass = self.initial_mass
            self.attitude = np.cross(np.cross(self.initial_position, self.initial_velocity), self.initial_position)
            self.attitude /= np.linalg.norm(self.attitude)

        def derivative(self, state):
            position = state[0]
            velocity = state[1]
            mass = state[2]
            if self.burn:
                thrust = self.max_thrust
            else:
                thrust = 0.0
            position_derivative = velocity
            velocity_derivative = -self.gravitational_parameter * position / np.linalg.norm(position) ** 3
            velocity_derivative += self.attitude * thrust / mass
            mass_derivative = -thrust / self.exhaust_velocity
            return np.array([position_derivative, velocity_derivative, mass_derivative])

        def step(self, delta_time):
            state = np.array([self.position, self.velocity, self.mass])
            k1 = delta_time * self.derivative(state)
            k2 = delta_time * self.derivative(state + k1 / 2)
            k3 = delta_time * self.derivative(state + k2 / 2)
            k4 = delta_time * self.derivative(state + k3)
            state += 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            self.position = state[0]
            self.velocity = state[1]
            self.mass = state[2]
            self.time += delta_time

        def energy(self):
            return np.linalg.norm(self.velocity) ** 2 / 2 - self.gravitational_parameter / np.linalg.norm(self.position)

        def semi_major_axis(self):
            return -self.gravitational_parameter / 2 / self.energy()

        def periapsis(self):
            return self.semi_major_axis() * (1 - self.eccentricity())

        def eccentricity_vector(self):
            return np.cross(self.velocity, np.cross(self.position, self.velocity)) / self.gravitational_parameter - self.position / np.linalg.norm(self.position)

        def eccentricity(self):
            return np.linalg.norm(self.eccentricity_vector())

        def true_anomaly(self):
            e = self.eccentricity_vector()
            r = self.position
            v = self.velocity
            vr = np.dot(r, v) / np.linalg.norm(r)
            true_anomaly = np.arccos(np.dot(e, r) / (np.linalg.norm(e) * np.linalg.norm(r)))
            if vr < 0:
                true_anomaly = 2 * np.pi - true_anomaly
            return true_anomaly

        def mean_anomaly(self):
            true_anomaly = self.true_anomaly()
            eccentricity = self.eccentricity()
            eccentric_anomaly = np.arctan2(np.sqrt(1 - eccentricity ** 2) * np.sin(true_anomaly), eccentricity + np.cos(true_anomaly))
            mean_anomaly = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)
            return mean_anomaly

        def integrate_to_eccentricity(self, eccentricity, delta_time):
            while self.eccentricity() < eccentricity:
                self.step(delta_time)
            while abs(self.eccentricity() - eccentricity) > 1e-8:
                delta_time *= 0.5
                if self.eccentricity() < eccentricity:
                    self.step(delta_time)
                else:
                    self.step(-delta_time)

    def burn_wait(self, angle_burn):
        vessel_angle = self.vessel.orbit.true_anomaly + self.vessel.orbit.argument_of_periapsis
        while vessel_angle > 2 * np.pi:
            vessel_angle -= 2 * np.pi
        while angle_burn > 2 * np.pi:
            angle_burn -= 2 * np.pi
        if angle_burn < vessel_angle:
            delta_true_anomaly = 2 * np.pi + angle_burn - vessel_angle
        else:
            delta_true_anomaly = angle_burn - vessel_angle

        position = np.array(self.vessel.position(self.ref_frame))
        velocity = np.array(self.vessel.velocity(self.ref_frame))
        normal = np.cross(position, velocity)
        attitude = um.normalize(np.cross(normal, position))
        attitude = um.vector_axis_angle(attitude, normal, delta_true_anomaly)

        time_wait = delta_true_anomaly / (2 * np.pi) * self.vessel.orbit.period
        return attitude, time_wait

    def execute(self):
        current_phase_time = int(self.time_start + self.time_launch)

        while current_phase_time > self.conn.space_center.ut:
            self.ui.update_phase_current_time(int(current_phase_time - self.conn.space_center.ut))

        current_phase_time = int(self.conn.space_center.ut + self.launch_time_offset)

        while self.vessel.control.throttle > 0.0:
            self.ui.update_phase_current_time(int(current_phase_time - self.conn.space_center.ut))

        self.ui.update_phase_current('Low Kerbin Orbit')
        self.ui.update_phase_next('TLI Burn')

        self.ui.update_message('')
        self.ui.update_phase_current_time(0)

        self.conn.space_center.physics_warp_factor = 0

        self.conn.space_center.warp_to(self.time_start + self.time_start_burn - 60)

        self.ui.update_phase_current('TLI Burn')
        self.ui.update_phase_next('Extract MEM')

        attitude, time_wait = self.burn_wait(self.angle_burn)

        uk.burn(self.conn, self.ui, self.ref_frame, self.burn_delta_velocity * attitude, delay=time_wait, delta_velocity_offset=-4, stopping_time=(10.0, 0.5, 10.0), change_ui_phases=True, phase_current='TLI Burn', phase_next='Extract MEM')

        self.ui.update_phase_current('Extract MEM')
        self.ui.update_phase_next('MCC')
        self.ui.update_next_action('Extract MEM')
        self.ui.update_phase_current_time(0)
        time.sleep(4)
        self.conn.space_center.warp_to(self.conn.space_center.ut + 300)

        while True:
            if self.ui.button.clicked:
                self.extract_lem()
                self.ui.button.clicked = False
                break

    def extract_lem(self):
        self.vessel.control.rcs = True
        self.vessel.control.sas = True
        self.vessel.control.activate_next_stage()
        time.sleep(1)
        self.vessel.control.activate_next_stage()
        time.sleep(1)
        self.vessel.control.forward = 1.0
        time.sleep(1)
        self.vessel.control.forward = 0.0
        time.sleep(1)
        self.vessel.control.forward = -1.0
        time.sleep(1)
        self.vessel.control.forward = 0.0

        while True:
            for x in self.conn.space_center.vessels:
                d = np.array(x.position(self.ref_frame)) - np.array(self.vessel.position(self.ref_frame))
                if 0 < um.magnitude(d) < 50:
                    self.conn.space_center.target_docking_port = x.parts.docking_ports[0]
                    break
            if self.conn.space_center.target_docking_port is not None:
                break
        time.sleep(0.5)
        self.vessel.control.sas_mode = self.vessel.control.sas_mode.target

        while um.vector_angle(np.array([0, 1, 0]), np.array(self.conn.space_center.target_docking_port.position(self.vessel.reference_frame))) > 0.03:
            time.sleep(0.2)

        time.sleep(4)

        self.vessel.control.forward = 1.0
        time.sleep(0.5)
        self.vessel.control.forward = 0.0

        current_mass = self.vessel.mass
        while True:
            if self.vessel.mass > current_mass + 1000:
                break
            time.sleep(1)

        self.vessel.control.toggle_action_group(2)
        self.vessel.control.forward = -1.0
        time.sleep(0.5)
        self.vessel.control.forward = 0.0
        time.sleep(1)

