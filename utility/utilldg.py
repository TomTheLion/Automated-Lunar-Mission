import numpy as np
import utility.utilmath as um


# approach phase targeting routine
def approach_phase_targeting_routine(mission):
    tmf = mission.tapm - mission.tapf
    tif = mission.tapi - mission.tapf

    apfgx = eq2102(tmf, tif, mission)
    apfgz = eq2104(tmf, tif, mission)

    apfg = []

    for i in range(5):
        apfg.append([apfgx[i], 0, apfgz[i]])

    mission.aptg = np.matmul(um.state_transition_matrix(0, mission.tapf), apfg)
    mission.apig = np.matmul(um.state_transition_matrix(mission.tapi, 0), mission.aptg)


# breaking phase targeting routine
def breaking_phase_targeting_routine(mission):
    rbrfg = mission.apig[0]
    vbrfg = mission.apig[1]
    unfbrfg = np.array([np.cos(mission.pbrf), 0, -np.sin(mission.pbrf)])
    abrfg = np.array([mission.fbrf / mission.mbrf * unfbrfg[0] - mission.surface_gravity, 0, mission.fbrf / mission. mbrf * unfbrfg[2]])
    mdbrf = -mission.fbrf / mission.vexh
    jbrfg = np.array([mission.jbrfg[0], 0, -mission.kj * abrfg[2] * mdbrf / mission.mbrf])
    sbrfg = np.array([mission.sbrfg[0], 0, mission.sbrfg[2]])

    mission.brfg = np.array([rbrfg, vbrfg, abrfg, jbrfg, sbrfg])
    mission.brtg = np.matmul(um.state_transition_matrix(0, mission.tbrf), mission.brfg)

    if mission.brig is None:
        mission.brig = np.matmul(um.state_transition_matrix(mission.tbri, 0), mission.brtg)


# equation 21.1
def eq2102(tmf, tif, mission):
    matrix_1 = np.array([[tmf ** 2 / 2, tmf ** 3 / 6, tmf ** 4 / 24], [tmf, tmf ** 2 / 2, tmf ** 3 / 6], [tif ** 2 / 2, tif ** 3 / 6, tif ** 4 / 24]])
    matrix_1 = np.linalg.inv(matrix_1)
    matrix_2 = np.array([[-1, -tmf, 1, 0, 0], [0, -1, 0, 1, 0], [-1, -tif, 0, 0, 1]])
    matrix_3 = np.array([mission.rapfgx, mission.vapfgx, mission.rapmg[0], mission.vapmg[0], mission.rapig[0]])
    matrix_4 = np.matmul(np.matmul(matrix_1, matrix_2), matrix_3)
    apfgx = np.array([mission.rapfgx, mission.vapfgx, matrix_4[0], matrix_4[1], matrix_4[2]])
    return apfgx


# equations 21.4 - 21.6
def eq2104(tmf, tif, mission):
    matrix_1 = np.array([[mission.tau ** 2 - mission.tau * tmf + tmf ** 2 / 2, tmf ** 3 / 6, tmf ** 4 / 24], [-mission.tau + tmf, tmf ** 2 / 2, tmf ** 3 / 6], [mission.tau ** 2 - mission.tau * tif + tif ** 2 / 2, tif ** 3 / 6, tif ** 4 / 24]])
    matrix_1 = np.linalg.inv(matrix_1)
    matrix_2 = np.array([mission.rapmg[2], mission.vapmg[2], mission.rapig[2]])
    matrix_3 = np.matmul(matrix_1, matrix_2)
    apfgz = np.array([matrix_3[0] * mission.tau ** 2, -matrix_3[0] * mission.tau, matrix_3[0], matrix_3[1], matrix_3[2]])
    return apfgz


# equations 21.9
def eq2109(t, matrix_2):
    matrix_1 = np.array([[-24 / t ** 3, -18 / t ** 2, - 6 / t, 24 / t ** 3, -6 / t ** 2], [72 / t ** 4, 48 / t ** 3, 12 / t ** 2, -72 / t ** 4, 24 / t ** 3]])
    matrix_3 = np.matmul(matrix_1, matrix_2)
    return matrix_3[0], matrix_3[1]


# p63, p64 guidance algorithm
def guide(tg, dt, r, v, l, mission):
    cgx = um.normalize(l)
    cgy = um.normalize(np.cross(l, r - l - mission.k * (v - np.cross(mission.body_rotation, r)) * tg / 4))
    cgz = um.normalize(np.cross(cgx, cgy))

    cg = np.array([cgx, cgy, cgz])
    rg = np.matmul(cg, r - l)
    vg = np.matmul(cg, v - np.cross(mission.body_rotation, r))

    a = mission.target[3][2]
    b = 6 * mission.target[2][2]
    c = 18 * mission.target[1][2] + 6 * vg[2]
    d = 24 * (mission.target[0][2] - rg[2])
    tg += dt

    while True:
        dt = -(a * tg ** 3 + b * tg ** 2 + c * tg + d) / (3 * a * tg ** 2 + 2 * b * tg + c)
        if abs(dt) > 60:
            dt *= 60 / abs(dt)
        tg = tg + dt
        if abs(dt) < 1 / 128:
            break

    tp = tg + mission.leadtime
    a = (3 * (tp / tg) ** 2 - 2 * (tp / tg)) * 12 / tg ** 2
    b = (4 * (tp / tg) ** 2 - 3 * (tp / tg)) * 6 / tg
    c = (2 * (tp / tg) ** 2 - (tp / tg)) * 6 / tg
    d = (6 * (tp / tg) ** 2 - 6 * (tp / tg) + 1)
    acg1 = a * (mission.target[0] - rg)
    acg2 = b * mission.target[1]
    acg3 = c * vg
    acg4 = d * mission.target[2]

    acg = acg1 + acg2 + acg3 + acg4
    g = -mission.gravitational_parameter / um.magnitude(r) ** 3 * r
    afcp = np.matmul(np.transpose(cg), acg) - g

    return afcp, tg, rg, vg


def landing_site(d, lat, lng):

    li = np.array([1, 0, 0])
    li = um.vector_axis_angle(li, np.array([0, 0, 1]), lng)
    li = d * um.vector_axis_angle(li, np.cross(li, np.array([0, 0, 1])), lat)

    return li

# p63 ignition algorithm
def p63_ignition_algorithm(mission):
    t = mission.time_initial
    tg = mission.guide_time_initial
    li = mission.landing_position_initial
    unfcp = -um.normalize(mission.velocity_initial)
    mission.target = mission.brtg

    while True:
        r, v = um.cse_rv(mission.position_initial, mission.velocity_initial, t, mission.gravitational_parameter)
        lt = um.vector_axis_angle(li, np.array([0, 1, 0]), mission.body_rotation[1] * t)
        for i in range(3):
            v += mission.aftrim * unfcp * mission.ttrim
            afcp, tg, rg, vg = guide(tg, 0, r, v, lt, mission)
            unfcp = um.normalize(afcp)
        dt = np.dot(np.array([rg[0] - mission.brig[0][0], rg[1] ** 2, rg[2] - mission.brig[0][2], um.magnitude(vg) - um.magnitude(mission.brig[1])]), mission.kgt) / (vg[2] + mission.kgt[0] * vg[0])

        if abs(dt) > 1:
            dt *= 1 / abs(dt)
        t -= dt
        if abs(dt) < 1 / 128:
            break

    t -= (mission.tullage + mission.ttrim)
    mission.time = t
    mission.guide_time = tg


# fourth order runge kutta
def runge_kutta(r0, v0, a, h, gm):
    def f(y):
        r = y[0]
        v = y[1]
        g = -gm / um.magnitude(r) ** 3 * r
        return np.array([v, a + g])

    y0 = np.array([r0, v0])

    k1 = h * f(y0)
    k2 = h * f(y0 + k1 * h / 2)
    k3 = h * f(y0 + k2 * h / 2)
    k4 = h * f(y0 + k3 * h)

    y1 = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return y1[0], y1[1]


def simulate(mission, target, time_limit, delta_time, afcp_limit, guide_time_flag):
    mission.target = target
    t0 = mission.time
    rg_out = []
    thrott_out = []
    while True:
        lt = um.vector_axis_angle(mission.landing_position_initial, np.array([0, 1, 0]), mission.body_rotation[1] * mission.time)
        afcp, mission.guide_time, rg, vg = guide(mission.guide_time, delta_time, mission.position, mission.velocity, lt, mission)
        maxafcp = mission.maxthrust / mission.mass
        if afcp_limit == -1:
            if um.magnitude(afcp) > maxafcp * mission.minthrottle:
                afcp = um.normalize(afcp) * maxafcp
            elif not mission.throttflag:
                mission.throttflag = True
                mission.throtta = mission.guide_time
        else:
            afcp = um.normalize(afcp) * afcp_limit
        mission.position, mission.velocity = runge_kutta(mission.position, mission.velocity, afcp, delta_time, mission.gravitational_parameter)
        mission.mass -= um.magnitude(afcp) * mission.mass / mission.vexh * delta_time
        mission.time += delta_time

        rg_out.append(rg)
        thrott_out.append([mission.guide_time, um.magnitude(afcp) / maxafcp])
        if guide_time_flag and mission.guide_time > time_limit:
            break
        if not guide_time_flag and mission.time - t0 > time_limit:
            break
    return rg_out, thrott_out


def simulation(mission):
    # set apig, aptg
    approach_phase_targeting_routine(mission)
    # set brig, brtg
    breaking_phase_targeting_routine(mission)
    # # set time, guide_time

    p63_ignition_algorithm(mission)
    # set position and velocity at ignition
    mission.time += mission.tullage
    mission.position, mission.velocity = um.cse_rv(mission.position_initial, mission.velocity_initial, mission.time, mission.gravitational_parameter)
    mission.mass = mission.mass_initial
    # # step forward, delta time is trim time

    rg_trim, thrott_trim = simulate(mission, mission.brtg, mission.ttrim, 1, mission.aftrim, False)
    rg_br, thrott_br = simulate(mission, mission.brtg, mission.tbrf, 1, -1, True)
    rg_ap, thrott_ap = simulate(mission, mission.aptg, mission.tapf, 0.2, -1, True)


def state_from_orbit(ra, rp, lat, lng, node_dir, vel_dir, gm):
    inc = lat + (np.pi / 2 - lat) * 0.05

    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)

    dt = 2 * np.pi * np.sqrt(a ** 3 / gm) / 6
    dt = 616.0809208241426

    r = np.array([1, 0, 0])
    r = um.vector_axis_angle(r, np.array([0, 0, 1]), lng)
    r = rp * um.vector_axis_angle(r, np.cross(r, np.array([0, 0, 1])), lat)

    longitude_of_node = lng + node_dir * abs(um.safe_asin(np.tan(lat) / np.tan(inc)))

    print('LAN', np.degrees(longitude_of_node))
    print('a', a)
    print('e', e)
    print('dt', dt)

    n = np.array([1, 0, 0])
    n = um.vector_axis_angle(n, np.array([0, 0, 1]), longitude_of_node)
    normal = vel_dir * um.normalize(np.cross(n, r))

    v = np.sqrt((1 + e) / (1 - e) * gm / a) * um.normalize(np.cross(normal, r))

    r, v = um.cse_rv(r, v, -dt, gm)

    return r, v
