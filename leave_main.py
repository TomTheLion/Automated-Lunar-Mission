import krpc
import numpy as np
import time
from classes.Interface import Interface
from classes.Mission import Mission
from classes.TransLunarInjection import TransLunarInjection
from classes.MidCourseCorrection import MidCourseCorrection
from classes.LunarOrbitInsertion import LunarOrbitInsertion
from classes.LunarDescentGuidance import LunarDescentGuidance

conn = krpc.connect(name='automated_lunar_mission')
planet = 'Kerbin'
moon = 'Mun'
ui = Interface(conn)

mission = Mission(conn, planet, moon, ui)
ui.update_phase_current('Solve Trajectory')
ui.update_phase_next('Launch to Orbit')

trans_lunar_injection = TransLunarInjection(mission)
trans_lunar_injection.solve_trajectory()
if trans_lunar_injection.convergence:
    ui.update_phase_current('Launch to Orbit')
    ui.update_phase_next('Low Kerbin Orbit')
    trans_lunar_injection.execute()

    ui.update_phase_current('MCC')
    ui.update_phase_next('Mun Arrival')
    ui.update_next_action('MCC')
    time.sleep(4)
    conn.space_center.warp_to(conn.space_center.ut + 30)

    mid_course_correction = MidCourseCorrection(mission)

    while True:
        if ui.button.clicked:
            mid_course_correction.calculate_correction(tli=True)
            mid_course_correction.execute_correction_rcs()
            ui.update_next_action('N/A')
            ui.button.clicked = False
            break

    ui.update_phase_current('Arrive at Mun')
    ui.update_phase_next('LOI Burn Start')
    ui.update_phase_current_time(0)
    time.sleep(4)

    conn.space_center.rails_warp_factor = 5
    while np.linalg.norm(conn.space_center.active_vessel.position(conn.space_center.bodies[moon].reference_frame)) > 2e6:
        time.sleep(0.1)

    conn.space_center.rails_warp_factor = 0

    lunar_orbit_insertion = LunarOrbitInsertion(mission)

    ui.update_phase_current('LOI Burn')
    ui.update_phase_next('CI Burn')
    lunar_orbit_insertion.circularize()

    ui.update_phase_current('CI Burn')
    ui.update_phase_next('Circ Burn')
    lunar_orbit_insertion.change_inclination()

    ui.update_phase_current('Circ Burn')
    ui.update_phase_next('Undock MEM')
    lunar_orbit_insertion.circularize()

    ui.update_phase_current('Undock MEM')
    ui.update_phase_next('LDO Burn')
    ui.update_next_action('Undock MEM')

    while True:
        if ui.button.clicked:
            conn.space_center.active_vessel.parts.docking_ports[0].undock()
            ui.update_next_action('LDO Burn')
            ui.button.clicked = False
            break

    ui.update_phase_current('LDO Burn')
    ui.update_phase_next('P63')

    while True:
        if ui.button.clicked:
            for engine in conn.space_center.active_vessel.parts.engines:
                if engine.part.tag == 'LDO Engine':
                    engine.active = True
            lunar_orbit_insertion.descent_orbit()
            ui.update_next_action('N/A')
            ui.button.clicked = False
            break

    conn.space_center.warp_to(conn.space_center.ut + conn.space_center.active_vessel.orbit.period * 0.375)

    lunar_descent_guidance = LunarDescentGuidance(mission)

    lunar_descent_guidance.execute()
