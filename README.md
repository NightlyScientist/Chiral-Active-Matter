---
authors: ["JImmy Gonzalez Nunez", "Subhaya Bose", "Arya Rajan"]
university: "University of California, Merced"
date: May 14, 2021
---

# Chrial Active Matter

Active spinner systems, a non-equilibrium system, exhibit unconventional
properties such as edge currents, frustration-induced melting (under pressure),
and periodicity in spatial ordering. [Von Zuiden et al][active-dimers] examined
the spatial characteristics in a system composed of representative diatomic
molecules, which have no active linear force and have an internal active torque.
The competition between the active rotation and short-range repulsive forces
leads to various distinct phases that are controlled by local packing fraction.
In this proposed work, we extend the work by considering a system of stiff rods.
We explore the effects of particle geometry on the emergent phase behaviour.

## Running the Code

### Requirements

- Python 3.8+
- numpy 1.20.2+
- numba 0.53.1+
- celluloid 0.2.0+
- matplotlib 3.4.2+
- argparse 3.2+
- numexpr 2.7.3+

### Run Python Script

The main file is called 'main.py'. Run the code by navigating to the _/src_
directory. Within the _/src_ directory, use

> Python main.py --help

to see all available command line flags. As an example, let's study the
behaviour of a system of $64$ rods at a packing fraction near $0.1$. We shall
_animate_ the system, calculate the _spatial correlation_ function, and tell our
script to print lots of information to the terminal (std):

> python main.py --phi 0.1 --verbose --animate --spatialCorr --frames 10 --spf
> 100 --color_orientation

here we are capturing only $10$ frames, with $100$ time steps between each
frame, and coloring the rods by their nematic orientation. We can save the
results (and supress showing animation at the end of the run) by calling the
_--save_ flag.

To change the save directory, use

> --dirname "path_to_your_directory"


[active-dimers]: https://doi.org/10.1073/pnas.1609572113 "Spatiotemporal order and emergent edge currents in active spinner materials"
