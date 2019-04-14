import numpy as np
from math import factorial


# conic state exrapolation, inputs: initial position, initial velocity, outputs: final position, final velocity
def cse_rv(r0, v0, dt, gm):
    rscale = np.linalg.norm(r0)
    vscale = np.sqrt(gm / rscale)
    r0s = r0 / rscale
    v0s = v0 / vscale
    dts = dt * vscale / rscale
    v2s = np.linalg.norm(v0) ** 2 * rscale / gm
    alpha = 2 - v2s
    armd1 = v2s - 1
    rvr0s = np.dot(r0, v0) / np.sqrt(gm * rscale)

    x = 0
    ratio = 1
    x2 = x * x
    z = alpha * x2
    s, c = sc_functions(z)
    x2c = x2 * c

    while abs(ratio) > 1e-10:
        f = x + rvr0s * x2c + armd1 * x * x2 * s - dts
        df = x * rvr0s * (1 - z * s) + armd1 * x2c + 1
        ratio = f / df
        x -= ratio
        x2 = x * x
        z = alpha * x2
        s, c = sc_functions(z)
        x2c = x2 * c

    lf = 1 - x2c
    lg = dts - x2 * x * s

    r1 = lf * r0s + lg * v0s
    ir1 = 1 / np.linalg.norm(r1)
    lfdot = ir1 * x * (z * s - 1)
    lgdot = 1 - x2c * ir1

    v1 = lfdot * r0s + lgdot * v0s

    return r1 * rscale, v1 * vscale


# conic state extrapolation, inputs: initial position, final position, outputs: initial velocity, final velocity
def cse_rr(r0, r1, dt, d, gm):
    rscale = np.linalg.norm(r0) + np.linalg.norm(r1)
    vscale = np.sqrt(gm / rscale)
    r0s = r0 / rscale
    r1s = r1 / rscale
    dts = dt * vscale / rscale

    r0m = np.linalg.norm(r0s)
    r1m = np.linalg.norm(r1s)

    theta = np.arccos(np.dot(r0s, r1s) / r0m / r1m)
    if d * np.cross(r0s, r1s)[1] > 0:
        theta = 2 * np.pi - theta

    a = np.sin(theta) * np.sqrt(r0m * r1m / (1 - np.cos(theta)))
    z = 0
    dz = 0
    while True:
        s, c = sc_functions(z)
        y = 1 + a * (z * s - 1) / np.sqrt(c)
        if y < 0:
            z -= dz / 2
            dz /= 2
        else:
            f = (y / c) ** 1.5 * s + a * np.sqrt(y) - dts
            if abs(z) < 1e-4:
                df = np.sqrt(2) / 40 * y ** 1.5 + a / 8 * (np.sqrt(y) + a * np.sqrt(1 / (2 * y)))
            else:
                df = (y / c) ** 1.5 * (1 / 2 / z * (c - 3 / 2 * s / c) + 3 / 4 * s * s / c) + a / 8 * (3 * s / c * np.sqrt(y) + a * np.sqrt(c / y))

            dz = -f / df

            if z + dz > 4 * np.pi ** 2:
                z = (z + 4 * np.pi ** 2) / 2
            else:
                z += dz

        if abs(f) < 1e-10:
            break

    lf = 1 - y / r0m
    lg = a * np.sqrt(y)
    lgdot = 1 - y / r1m
    v0 = 1 / lg * (r1s - lf * r0s)
    v1 = 1 / lg * (lgdot * r1s - r0s)

    return v0 * vscale, v1 * vscale


# calculate magnitude of vector
def magnitude(vector):
    return np.linalg.norm(vector)


# normalize vector
def normalize(vector):
    return vector / np.linalg.norm(vector)


# safe asin handles situations where x should be equal to +/- 1 but due to numerical error exceeds those values
def safe_asin(x):
    if x >= 1:
        return np.pi / 2
    if x <= -1:
        return -np.pi / 2
    return np.arcsin(x)


# safe acos handles situations where x should be equal to +/- 1 but due to numerical error exceeds those values
def safe_acos(x):
    if x >= 1:
        return 0
    if x <= -1:
        return np.pi
    return np.arccos(x)


# stumpff s and c functions
def sc_functions(z):
    az = abs(z)
    if az < 1e-4:
        s = (1 - z * (0.05 - z / 840)) / 6
        c = 0.5 - z * (1 - z / 30) / 24
    else:
        saz = np.sqrt(az)
        if z > 0:
            x = saz
            s = (saz - np.sin(x)) / (saz * az)
            c = (1 - np.cos(x)) / az
        else:
            x = np.exp(saz)
            s = (0.5 * (x - 1 / x) - saz) / (saz * az)
            c = (0.5 * (x + 1 / x) - 1) / az
    return s, c


# state transition matrix
def state_transition_matrix(t1, t0):
    dt = t1 - t0
    xlist = []
    ylist = []
    for i in range(5):
        for j in range(5):
            if j >= i:
                yval = j - i
                xval = dt ** yval / factorial(yval)
            else:
                xval = 0
            xlist.append(xval)
        ylist.append(xlist)
        xlist = []
    return np.vstack(ylist)


# angle between two vectors
def vector_angle(a, b):
    return safe_acos(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


# vector rotation
def vector_axis_angle(vector, axis, angle):
    axis = axis / np.linalg.norm(axis)
    return vector * np.cos(angle) + np.cross(axis, vector) * np.sin(angle) + axis * np.dot(axis, vector) * (1 - np.cos(angle))


# projection of a vector onto the plane normal to another
def vector_exclude(normal_vector, vector):
    return vector - np.dot(vector, normal_vector) / np.dot(normal_vector, normal_vector) * normal_vector


# convert to/from ksp coordinate system
def convert_yz(arr):
    return np.array([arr[0], arr[2], arr[1]])
