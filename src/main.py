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
parser.add_argument("--phi", type=float, default=0.1, help="particle packing fraction")
parser.add_argument(
    "--number_particles", type=int, default=64, help="number of rods to simulate"
)
parser.add_argument(
    "--animate",
    action="store_true",
    help="if flag is called, program will animate system",
)
# parser.add_argument("--show", action="store_true", help="if flag is called, program will display figures and animation (calls plt.show())")
parser.add_argument(
    "--save",
    action="store_true",
    help="save animation and any figures that are created",
)
parser.add_argument("--fps", type=int, default=7, help="fps of saved animation")
parser.add_argument(
    "--spf",
    type=int,
    default=100,
    help="specify the number of time steps between captured frames",
)
parser.add_argument(
    "--frames", type=int, default=20, help="specify the total number of captured frames"
)
parser.add_argument(
    "--spatialCorr",
    action="store_true",
    help="calculate the spatial distribution funciton at the end of the simulation",
)
parser.add_argument(
    "--color_orientation",
    action="store_true",
    help="if used, will using orientation as color in animation instead of instantaneous torque",
)
parser.add_argument("--dirname", default="./", help="save location of animation")
parser.add_argument(
    "--initMethod",
    default="random",
    help="initialization method ['random','square', 'hexaxtic', 'two_particles']",
)
parser.add_argument(
    "--verbose", action="store_true", help="show messages, warnings, etc"
)
parser.add_argument(
    "--lightTheme",
    action="store_true",
    help="uses a light theme instead of a dark theme",
)
parser.add_argument(
    "--CoM", action="store_true", help="plots only the center of mass of each rod"
)
args = parser.parse_args()


# number of particles
nparticles = args.number_particles

# packing fraction
packFrac = args.phi

# time step (using 'forward Euler' integration)
dt = 0.002

# non-dimensionalized intrinsic angular speed
alpha = 0.05

# diameter of beads
sigma = 1.0

# friction coefficient
rat_gam = 1.0 * sigma ** 2

# rotational and translational diffusion constants
D_T = 0.05
D_R = 0.05
D_T = np.sqrt(2.0 * D_T * dt)
D_R = np.sqrt(2.0 * D_R * dt)

# interaction cut off dist
cut_off = 5 * sigma

# inverse screening length
lam = 1 / sigma

# marker size in animation
psize = 20 * np.pi * sigma

# box size
BoxL = int(np.sqrt(np.pi * nparticles * sigma ** 2 / packFrac))

# steps per frame
spf = args.spf

# number of frames
frames = args.frames

# initial particles
ptcls = np.array([])
if args.initMethod.lower() == "random":
    ptcls = np.array(initial_conditions.random_initialize(nparticles, BoxL, sigma))
elif args.initMethod.lower() == "two_particles":
    ptcls = np.array(
        initial_conditions.two_paricle_initialize(
            BoxL / 2.0,
            BoxL / 2.0,
            BoxL / 2.0 + 1.5,
            BoxL / 2.0,
            np.pi / 2,
            np.pi / 2,
            BoxL,
            sigma,
        )
    )
elif args.initMethod.lower() == "square":
    ptcls = np.array(initial_conditions.square_initialize(nparticles, BoxL, sigma))
elif args.initMethod.lower() == "hexatic":
    if args.verbose:
        print("  ! using experimental initialization method\n")
    ptcls = np.array(initial_conditions.hexatic_initialize(nparticles, BoxL, sigma))

if args.verbose:
    print("  ! system parameters:")
    print("  ------------------")
    print(
        f"   > actual packing fraction: {metrics.packingFraction(BoxL, nparticles, sigma):.3f}\n"
    )
    print(f"   > number of particles: {nparticles}\n")
    print(f"   > box size: {BoxL}\n")
    print(f"   > frames, spf, and fps: {args.frames}, {args.spf}, {args.fps}\n")
    print(f"   > using {args.initMethod.lower()} initialization\n")
    if args.save and args.animate:
        print(f"   > saving animation to {args.dirname}\n")
    elif args.animate and not args.save:
        print(f"   > animation will be made and displayed, but not saved\n")
    print("  ------------------\n")

idxs = np.array(list(range(len(ptcls))))

# cho0se theme
cmap_orientation = "jet"
cmap_torque = "seismic"
if args.lightTheme:
    plt.style.use("Solarize_Light2")
    cmap_orientation = "winter"
    cmap_torque = "brg"
else:
    plt.style.use("dark_background")

# create figure and axis
fig, ax = plt.subplots()

ax.axis("off")
ax.set_xlim(0, int(BoxL))
ax.set_ylim(0, int(BoxL))

camera = Camera(fig)

nframes = 1

# model step (main loop)
for t in range(spf * frames):
    # interact pari-wise
    for a, b in interacting_pairs(ptcls, idxs, cut_off, BoxL):
        interact(ptcls[a], ptcls[b], lam, BoxL, dt)

    # update array of rods
    ptcls = np.array(
        [update(ptcl, sigma, alpha, dt, BoxL, D_T, D_R, rat_gam) for ptcl in ptcls]
    )

    # take snapshots for animation
    if np.mod(t, spf) == 0:
        if args.CoM:
            if args.color_orientation:
                ax.scatter(
                    [ptcl.cmx for ptcl in ptcls],
                    [ptcl.cmy for ptcl in ptcls],
                    s=psize,
                    c=[np.mod(ptcl.theta, np.pi) for ptcl in ptcls],
                    vmin=0,
                    vmax=np.pi,
                    cmap=cmap_orientation,
                    marker="o",
                )
            else:
                ax.scatter(
                    [ptcl.cmx for ptcl in ptcls],
                    [ptcl.cmy for ptcl in ptcls],
                    s=psize,
                    c=[ptcl.omega for ptcl in ptcls],
                    vmin=-2,
                    vmax=2,
                    cmap=cmap_torque,
                    marker="o",
                )
        else:
            for ptcl in ptcls:
                if args.color_orientation:
                    ax.scatter(
                        ptcl.rsx,
                        ptcl.rsy,
                        s=psize,
                        c=[np.mod(ptcl.theta, np.pi)] * 4,
                        vmin=0,
                        vmax=np.pi,
                        cmap=cmap_orientation,
                        marker="o",
                    )
                else:
                    ax.scatter(
                        ptcl.rsx,
                        ptcl.rsy,
                        s=psize,
                        c=[ptcl.omega] * 4,
                        vmin=-2,
                        vmax=2,
                        cmap=cmap_torque,
                        marker="o",
                    )

        camera.snap()
        if args.verbose:
            print(
                "  >>>  creating frame ",
                nframes,
                " of ",
                args.frames,
                " frames",
                end="\r",
            )
        nframes += 1


# with open("data_packfuck_phi.data", "a") as f:
#     f.write(
#         f"{metrics.packingFraction(BoxL, nparticles, sigma)},{np.mean([alpha + p.omega for p in ptcls])}\n"
#     )

# create animation
plt.title(rf"$\phi = {metrics.packingFraction(BoxL, nparticles, sigma):.2f}$")

if args.animate:
    anim = camera.animate()
    if args.save:
        anim.save(
            f"{args.dirname}/animation_nparticles_{nparticles}_omega{alpha:.2f}_phi_{packFrac:.3f}.gif",
            writer="ffmpeg",
            fps=args.fps,
        )
    else:
        plt.show()


if args.spatialCorr:
    grfig, grax = plt.subplots()
    gravg, radii, _ = metrics.TwoPointCorr(
        [p.cmx for p in ptcls],
        [p.cmy for p in ptcls],
        BoxL,
        max([BoxL / 2.8, 5 * sigma]),
        sigma,
    )
    grax.plot(radii, gravg, "g+")
    grax.set_ylabel(r"$g(r)$")
    grax.set_xlabel("radius")
    if args.save:
        grfig.savefig(f"{args.dirname}/spatial-correlation.png")
    else:
        plt.show()
