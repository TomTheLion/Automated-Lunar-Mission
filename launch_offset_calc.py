import krpc

conn = krpc.connect(name='tli-burn')
moon = conn.space_center.bodies['Mun']
planet = conn.space_center.bodies['Kerbin']
vessel = conn.space_center.active_vessel
ref_frame = planet.non_rotating_reference_frame

throttle = conn.get_call(getattr, vessel.control, 'throttle')

expr = conn.krpc.Expression.less_than(conn.krpc.Expression.call(throttle), conn.krpc.Expression.constant_float(1))
event = conn.krpc.add_event(expr)

with conn.stream(getattr, vessel, 'thrust') as thrust:
    with thrust.condition:
        while not thrust():
            thrust.wait()

print('engine start')

start_time = conn.space_center.ut

with event.condition:
    event.wait()
    print('engine cuttoff')

end_time = conn.space_center.ut
delta_time = end_time - start_time

angle_from_longitude_of_ascending_node = vessel.orbit.argument_of_periapsis + vessel.orbit.true_anomaly

print("delta time: ", delta_time)
print("angle from LAN: ", angle_from_longitude_of_ascending_node)
print("vessel mass: ", vessel.mass)
