import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
from particle import interact, update, interacting_pairs
import initial_conditions
import metrics
import time
import argparse
from metrics import cut_arrow

parser = argparse.ArgumentParser()
parser.add_argument("--phi", type=float, default=0.1)
args = parser.parse_args()


packFrac = args.phi

dt = 0.002
alpha = 0.05
sigma = 1.0
rat_gam = 1.0 * sigma ** 2
D_T = 0.05
D_R = 0.05

D_T = np.sqrt(2.0 * D_T * dt)
D_R = np.sqrt(2.0 * D_R * dt)

# print(f"D_T = {D_T:.2f}")
# print(f"D_R = {D_R:.2f}")

twopi = np.pi * 2

cut_off = 4 * sigma  # interaction cut off dist
lam = 1 / sigma  # inverse screening length
psize = 20 * np.pi * sigma  # marker size in animation

# number of particles
nparticles = 81

# box size
BoxL = int(np.sqrt(np.pi * nparticles * sigma ** 2 / packFrac))

# steps per frame
spf = 100

# number of frames
frames = 50

# initial particles
ptcls = np.array(initial_conditions.random_initialize(nparticles, BoxL, sigma))

# ptcls = np.array(initial_conditions.two_paricle_initialize(BoxL/2., BoxL/2., BoxL/2.+ 1.5, BoxL/2., np.pi/2, np.pi/2, BoxL, sigma))

# ptcls = np.array(initial_conditions.threesome_paricle_initialize(BoxL/2., BoxL/2., BoxL/2.+ 1.5, BoxL/2., BoxL/2., BoxL/2. + 2.5, np.pi/2, np.pi/2, 0., BoxL, sigma))

idxs = np.array(list(range(len(ptcls))))

# create figure for animation
# plt.style.use("Solarize_Light2")
plt.style.use("dark_background")
fig, ax = plt.subplots()
ax.axis("off")
ax.scatter([], [], s=psize, vmin=0, vmax=twopi)
cx, cy, vx, vy = zip(*((p.cmx, p.cmy, p.vx, p.vy) for p in ptcls))
# ax.quiver(cx, cy, vx, vy, scale=20.)
ax.set_xlim(0, int(BoxL))
ax.set_ylim(0, int(BoxL))
camera = Camera(fig)

startt = time.time()

nframes = 1

# model step (main loop)
for t in range(spf * frames):
    for a, b in interacting_pairs(ptcls, idxs, cut_off, BoxL):
        interact(ptcls[a], ptcls[b], lam, BoxL, dt)

    ptcls = np.array(
        [update(ptcl, sigma, alpha, dt, BoxL, D_T, D_R, rat_gam) for ptcl in ptcls]
    )

    if np.mod(t, spf) == 0:
        for ptcl in ptcls:
            # ax.scatter(
            #     ptcl.rsx,
            #     ptcl.rsy,
            #     s=psize,
            #     c=[np.mod(ptcl.theta, np.pi)] * 4,
            #     vmin=0,
            #     vmax=twopi,
            #     cmap="tab10",
            #     marker="h"
            # )
            ax.scatter(
                ptcl.rsx,
                ptcl.rsy,
                s=psize,
                c=[ptcl.omega] * 4,
                vmin=-2,
                vmax=2,
                cmap="seismic",
                marker="h",
            )

        # ax.scatter(
        # [ptcl.cmx for ptcl in ptcls],
        # [ptcl.cmy for ptcl in ptcls],
        # s=psize,
        # c=[(alpha + ptcl.omega) for ptcl in ptcls],
        # vmin=-2,
        # vmax=2,
        # cmap="seismic",
        # marker="h",
        # )
        # cx, cy, vx, vy = zip(*((p.cmx, p.cmy, p.vx, p.vy) for p in ptcls))
        # ax.quiver(cx, cy, vx, vy, np.hypot(vx,vy), pivot"mid", units="xy", color="green")
        # scale=np.hypot(ptcl.vx, ptcl.vy))
        # ax.scatter(ptcl.cmx, ptcl.cmy, s=psize, c="black", marker="h")
        camera.snap()
        print("frames ", nframes, end="\r")
        nframes += 1


endt = time.time()
# print("time = ", endt - startt, " (s)")

# with open("data_packfuck_phi.data", "a") as f:
#     f.write(
#         f"{metrics.packingFraction(BoxL, nparticles, sigma)},{np.mean([alpha + p.omega for p in ptcls])}\n"
#     )

# create animation
plt.title(rf"$\phi = {metrics.packingFraction(BoxL, nparticles, sigma):.2f}$")
anim = camera.animate()
# anim.save(f"ensemble_color_omega_phi_{packFrac:.3f}.mp4", writer="ffmpeg", fps=7)
# print(metrics.packingFraction(BoxL, nparticles, sigma))


grfig, grax = plt.subplots()
gravg, radii, _ = metrics.TwoPointCorr(
    [p.cmx for p in ptcls], [p.cmy for p in ptcls], BoxL, 5, sigma
)
grax.plot(radii, gravg, "k+")
grax.set_ylabel(r"$g(r)$")
grax.set_xlabel("radius")

plt.show()
