import numpy as np
import numexpr as ne

def packingFraction(BoxL, nParticles, sigma):
    return np.pi*nParticles*sigma**2/BoxL**2

def phi_6():
    pass

def TwoPointCorr(x, y, S, rMax, dr):
    print(f"  ! calculating g(r) using S, rmax, dr of ({S}, {rMax:.4f}, {dr:.4f})\n")
    bools1 = ne.evaluate("x > rMax")
    bools2 = ne.evaluate("x < (S - rMax)")
    bools3 = ne.evaluate("y > rMax")
    bools4 = ne.evaluate("y < (S - rMax)")

    (interiorIndices,) = np.where(bools1 * bools2 * bools3 * bools4)
    nInterPart = len(interiorIndices)

    if nInterPart < 1:
        print("  !! rMax is too large!\n")
        return [], [], []


    edges = np.arange(0.0, rMax + 1.1 * dr, dr)
    nIncrements = len(edges) - 1
    g = np.zeros([nInterPart, nIncrements])
    radii = np.zeros(nIncrements)
    numberDensity = len(x) / S ** 2

    for p in range(nInterPart):
        index = interiorIndices[p]
        d = np.sqrt((x[index] - x) ** 2 + (y[index] - y) ** 2)
        d[index] = 2 * rMax

        (results, bins) = np.histogram(d, bins=edges, normed=False)
        g[p, :] = results / numberDensity

    gAvg = np.zeros(nIncrements)
    for i in range(nIncrements):
        radii[i] = (edges[i] + edges[i + 1]) / 2.0
        rOuter = edges[i + 1]
        rInner = edges[i]
        gAvg[i] = np.mean(g[:, i]) / (np.pi * (rOuter ** 2 - rInner ** 2))

    return (gAvg, radii, interiorIndices)

def cut_arrow(c, v, BoxL):
    if v + c >= BoxL:
        return BoxL - c
    elif v + c >= -BoxL:
        return -c
    else:
        return c+v