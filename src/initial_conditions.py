import particle
import numpy as np
from numpy.random import uniform

def two_paricle_initialize(x1, y1, x2, y2, o1, o2, BoxL, sigma):
    ptcls = []
    ptcls.append(particle.Particle(x1, y1, sigma,o1, BoxL))
    ptcls.append(particle.Particle(x2, y2, sigma,o2, BoxL))
    return ptcls

def threesome_paricle_initialize(x1, y1, x2, y2, x3, y3, o1, o2, o3, BoxL, sigma):
    ptcls = []
    ptcls.append(particle.Particle(x1, y1, sigma,o1, BoxL))
    ptcls.append(particle.Particle(x2, y2, sigma,o2, BoxL))
    ptcls.append(particle.Particle(x3, y3, sigma,o3, BoxL))
    return ptcls
def random_initialize(nPartcles, BoxL, sigma):
    ptcls = [
        particle.Particle(
            uniform() * BoxL, uniform() * BoxL, sigma, uniform() * 2 * np.pi, BoxL
        )
        for _ in range(nPartcles)
    ]
    return ptcls


def square_initialize(nparticles, BoxL, sigma):
    nrow = int(2 * np.sqrt(nparticles))
    ncolumn = nparticles // nrow
    # cellcount = 0
    ptcls = []
    #   xPos = np.zeros(nparticles)
    #   yPos = np.zeros(nparticles)
    #   xsigma = 0.0
    for P in range(ncolumn):
        for Q in range(nrow):
            xsigma = P * BoxL / ncolumn + ncolumn / 2.0
            ysigma = Q * BoxL / nrow + nrow / 2.0
            # print(f"(x,y) = ({xsigma}, {ysigma})")
            ptcls.append(particle.Particle(xsigma, ysigma, sigma, np.pi, BoxL))
    return ptcls


def hexatic_initialize(nparticles, BoxL, sigma):
    nrow = int(2 * np.sqrt(nparticles))
    ncolumn = nparticles // nrow
    ptcls = []
    for P in range(ncolumn):
        if P % 2 == 0:
            for Q in range(nrow):
                if Q % 2 == 0:
                    xsigma = P * BoxL / ncolumn + ncolumn / 2.0
                    ysigma = 2 * Q * BoxL / nrow + nrow / 2.0
                    ptcls.append(particle.Particle(xsigma, ysigma, sigma, np.pi, BoxL))
        else:
            for Q in range(nrow):
                if Q % 2 == 1:
                    xsigma = P * BoxL / ncolumn + ncolumn / 2.0
                    ysigma = 2 * Q * BoxL / nrow + nrow / 2.0
                    ptcls.append(particle.Particle(xsigma, ysigma, sigma, np.pi, BoxL))
    return ptcls
