import kwant
from matplotlib import pyplot
import numpy as np
a = 1
t = 1.0
L = 30
W = 10

lat = kwant.lattice.square(a)
syst = kwant.Builder(kwant.TranslationalSymmetry(lat.vec((0,W)))
syst[(lat(x, y) for x in range(L) for y in range(W))] = 4 * t  # makes lattice
syst[lat.neighbors(1)] = -t
syst = kwant.wraparound.wraparound(syst)
