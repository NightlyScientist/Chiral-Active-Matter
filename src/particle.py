import numpy as np
from numpy import sin, cos, array, sqrt, exp, pi, fromiter
from numpy.random import normal
from itertools import combinations
from numba import njit, jit
from numba import types
from numba.extending import typeof_impl
import numexpr as ne

class Particle:
    def __init__(self, x, y, sigma, theta, BoxL, omega=1.0, vx=1., vy=1.):
        self.cmx = np.mod(x, BoxL)
        self.cmy = np.mod(y, BoxL)
        self.vx = vx
        self.vy = vy
        self.theta = theta
        self.omega = omega
        self.rsx = np.mod(array([-1.5, -0.5, 0.5, 1.5]) * sigma * cos(theta) + x, BoxL)
        self.rsy = np.mod(array([-1.5, -0.5, 0.5, 1.5]) * sigma * sin(theta) + y, BoxL)
        self.fx = array([0.0, 0.0, 0.0, 0.0])
        self.fy = array([0.0, 0.0, 0.0, 0.0])


class ParticleType(types.Type):
    def __init__(self):
        super(ParticleType, self).__init__(name="Particle")


config = {"nogil": True, "fastmath": True, "cache": True, "nopython": True}

particle_type = ParticleType()


@typeof_impl.register(Particle)
def typeof_index(val, c):
    return particle_type


def distance(ptcla, ptclb, cut_off, BoxL):
    # fixed dist to account for periodic walls
    return (metric(ptcla.cmx - ptclb.cmx, BoxL)) ** 2 + (
        metric(ptcla.cmy - ptclb.cmy, BoxL)
    ) ** 2 <= cut_off ** 2


def interacting_pairs(ptcls, idxs, cut_off, BoxL):
    # get unique pairs of idices
    cmbtns = combinations(idxs, 2)
    # remove pairs whose dist > cut_off
    return (x for x in cmbtns if distance(ptcls[x[0]], ptcls[x[1]], cut_off, BoxL))


    # return filter(lambda pt: (ptcls[pt[1]].cmx - ptcls[pt[1]].cmx)**2
    # + (ptcls[pt[0]].cmy - ptcls[pt[1]].cmy)**2 <= cut_off**2, cmbtns)


def unitvector(v):
    return v / (v ** 2).sum() ** 0.5


@jit(**config)
def crossprod_3(rx, ry, fx, fy):
    return rx * fy - ry * fx


def update(ptcl, sigma, alpha, dt, BoxL, D_T, D_R, rat_gam):
    vx = ptcl.fx.sum()
    vy = ptcl.fy.sum()
    cmMotionx = (ptcl.cmx + vx * dt) #+ D_T * normal(0, 1, 1)
    cmMotiony = (ptcl.cmy + vy * dt) #+ D_T * normal(0, 1, 1)
    torque = (
        rat_gam
        * crossprod_3(ptcl.rsx - ptcl.cmx, ptcl.rsy - ptcl.cmy, ptcl.fx, ptcl.fy).sum()
    )
    theta = (ptcl.theta + (torque + alpha) * dt) #+ D_R * normal(0, 1, 1)
    # cmMotionx = ptcl.cmx
    # cmMotiony = ptcl.cmy
    # print(ptcl.rsx)
    # print(ptcl.rsy)
    # print(ptcl.fx)
    # print(ptcl.fy)
    # print(torque, " <? ", alpha)
    return Particle(
        cmMotionx, cmMotiony, sigma, np.mod(theta, 2.0 * pi), BoxL, omega=torque, vx=vx, vy=vy
    )


def interact(ptcla, ptclb, lamb, BoxL, dt):
    nodeRepulsion(ptcla, ptclb, lamb, 0, BoxL, dt)
    nodeRepulsion(ptcla, ptclb, lamb, 1, BoxL, dt)
    nodeRepulsion(ptcla, ptclb, lamb, 2, BoxL, dt)
    nodeRepulsion(ptcla, ptclb, lamb, 3, BoxL, dt)
    # nodeRepulsion(ptcla, ptclb, lamb, 5)


# @njit
@jit(**config)
def metric(x, BoxL):
    if x > BoxL / 2.0:
        return x - BoxL
    elif x < -BoxL / 2.0:
        return x + BoxL
    else:
        return x


@jit(**config)
def force(r, lamb, dt):
    if r < 1.0:
        return 1 / dt * (1 - r)
    else:
        return exp(-lamb * r) * (lamb * r + 1) / (r ** 2)


# @njit
@jit(**config)
def separation(rx, ry):
    r0 = sqrt(rx ** 2 + ry ** 2)
    hatr0x = rx / r0
    hatr0y = ry / r0
    return r0, hatr0x, hatr0y

    # f0 = exp(-lamb * r0) * (lamb * r0 + 1) / (r0 ** 2)  # TODO :add b coeff
    # f0 = array([force(r, lamb, dt) for r in r0])

def nodeRepulsion(ptcla, ptclb, lamb, idx, BoxL, dt):
    rx = -(ptclb.rsx - ptcla.rsx[idx])
    ry = -(ptclb.rsy - ptcla.rsy[idx])
    rx = fromiter((metric(x, BoxL) for x in rx), rx.dtype, count=4)
    ry = fromiter((metric(y, BoxL) for y in ry), ry.dtype, count=4)
    r0, hatr0x, hatr0y = separation(rx, ry)
    f0 = fromiter((force(r, lamb, dt) for r in r0), r0.dtype, count=4)
    f0xs = ne.evaluate("f0 * hatr0x")
    f0ys = ne.evaluate("f0 * hatr0y")
    ptcla.fx[idx] += f0xs.sum()
    ptclb.fx -= f0xs
    ptcla.fy[idx] += f0ys.sum()
    ptclb.fy -= f0ys
