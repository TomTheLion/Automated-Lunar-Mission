import krpc
import numpy as np
import time
from classes.Interface import Interface
from classes.Mission import Mission
from classes.MidCourseCorrection import MidCourseCorrection
from classes.LunarReturn import LunarReturn
from classes.Rendezvous import Rendezvous


conn = krpc.connect(name='automated_lunar_mission')
planet = 'Kerbin'
moon = 'Mun'
ui = Interface(conn)

mission = Mission(conn, planet, moon, ui)

ui.update_phase_current('Launch to Orbit')
ui.update_phase_next('Rendezvous')

ui.update_next_action('N/A')

conn.space_center.active_vessel.parts.engines[0].active = False
conn.space_center.active_vessel.parts.engines[0].part.engine.thrust_limit = 0.40
conn.space_center.active_vessel.control.rcs = True
conn.space_center.active_vessel.parts.controlling.rcs.enabled = True

while conn.space_center.active_vessel.control.throttle < 1.0:
    pass

start_time = conn.space_center.ut

while conn.space_center.ut < start_time + 30:
    ui.update_phase_current_time(int(start_time + 30 - conn.space_center.ut))

conn.space_center.active_vessel.parts.engines[0].active = True

current_phase_time = conn.space_center.ut + 215
while conn.space_center.active_vessel.control.throttle > 0.0:
    ui.update_phase_current_time(int(current_phase_time - conn.space_center.ut))

ui.update_phase_current_time(0)
ui.update_phase_current('Rendezvous')
ui.update_phase_next('Undock MEM')

ui.update_next_action('N/A')

conn.space_center.active_vessel.parts.engines[0].part.engine.thrust_limit = 1.00
rend = Rendezvous(mission)
rend.execute_rendezvous()

ui.update_phase_current('Undock MEM')
ui.update_phase_next('TEI Burn')
ui.update_next_action('Undock MEM')

while True:
    if ui.button.clicked:
        conn.space_center.active_vessel.parts.docking_ports[0].undock()
        ui.update_next_action('N/A')
        ui.button.clicked = False
        break

ui.update_phase_current('TEI Burn')
ui.update_phase_next('MCC')
ui.update_next_action('MCC')

lr = LunarReturn(mission)
lr.calculate_return()
lr.execute_correction()

ui.update_phase_current('MCC')
ui.update_phase_next('Jettison SM')
ui.update_next_action('MCC')

mid_course_correction = MidCourseCorrection(mission)

while True:
    if ui.button.clicked:
        mid_course_correction.calculate_correction()
        mid_course_correction.execute_correction()
        ui.button.clicked = False
        break

print(conn.space_center.active_vessel.velocity(conn.space_center.bodies[planet].non_rotating_reference_frame))

while True:
    if ui.button.clicked:
        mid_course_correction.calculate_correction()
        mid_course_correction.execute_correction_rcs()
        ui.update_next_action('Jettison SM')
        ui.button.clicked = False
        break

ui.update_phase_current('Jettison SM')
ui.update_phase_next('Deploy Chute')

conn.space_center.warp_to(conn.space_center.ut + conn.space_center.active_vessel.orbit.time_to_periapsis - 120)

while True:
    if ui.button.clicked:
        conn.space_center.active_vessel.control.toggle_action_group(3)
        ui.update_next_action('Deploy Chute')
        ui.button.clicked = False
        break

conn.space_center.active_vessel.control.sas = True
time.sleep(2 )
conn.space_center.active_vessel.control.sas_mode = conn.space_center.active_vessel.control.sas_mode.retrograde

ui.update_phase_current('Deploy Chute')
ui.update_phase_next('Land')

while True:
    if ui.button.clicked:
        conn.space_center.active_vessel.control.toggle_action_group(4)
        ui.update_next_action('N/A')
        ui.button.clicked = False
        break

ui.update_phase_current('Land')
ui.update_phase_next('')

while conn.space_center.active_vessel.situation.name is not 'splashed' or 'landed':
    pass
