import numpy as np


def burn(conn, ui, ref_frame, delta_velocity, delay=0, delta_velocity_offset=0, stopping_time=(0.5, 0.5, 0.5), change_ui_phases=False, phase_current='', phase_next='', next_action='N/A', main_engine=True, specified_thrust=False, specified_isp=False, facing=1):
    vessel = conn.space_center.active_vessel

    burn_start = conn.space_center.ut + delay

    ft = vessel.max_thrust
    if specified_thrust:
        ft = specified_thrust

    isp = vessel.specific_impulse
    if specified_thrust:
        isp = specified_isp

    if delta_velocity is not False:
        dv = np.linalg.norm(delta_velocity) + delta_velocity_offset
        dm = vessel.mass * (1 - 1 / np.exp(dv / 9.80665 / isp))
        dt = dm / (ft / (9.80665 * isp))

    udv = delta_velocity / np.linalg.norm(delta_velocity) * facing

    vessel.auto_pilot.engage()
    vessel.auto_pilot.reference_frame = ref_frame
    vessel.auto_pilot.stopping_time = stopping_time
    vessel.auto_pilot.target_direction = udv

    settle_time = conn.space_center.ut + 5

    while np.linalg.norm(vessel.angular_velocity(ref_frame)) > 0.01 or conn.space_center.ut < settle_time or burn_start > conn.space_center.ut:
        if delay != 0:
            ui.update_phase_current_time(abs(int(burn_start - conn.space_center.ut)))

    burn_start = conn.space_center.ut
    if main_engine:
        vessel.control.throttle = 1.00 * facing
    else:
        vessel.control.rcs = True
        vessel.control.forward = 1.00 * facing

    if change_ui_phases:
        ui.update_phase_current(phase_current)
        ui.update_phase_next(phase_next)
        ui.update_next_action(next_action)

    while burn_start + dt > conn.space_center.ut:
        ui.update_phase_current_time(abs(int(burn_start + dt - conn.space_center.ut)))

    if main_engine:
        vessel.control.throttle = 0.00
    else:
        vessel.control.rcs = False
        vessel.control.forward = 0.00

    vessel.auto_pilot.disengage()
