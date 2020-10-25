import kwant
from matplotlib import pyplot
import numpy as np
a = 1
t = 1.0
L = 30
W = 10

def plot_conductance(syst, energies):
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()

lat = kwant.lattice.square(a)

syst = kwant.Builder(kwant.TranslationalSymmetry(lat.vec((0,L))))
syst[(lat(x, y) for x in range(L) for y in range(W))] = 4 * t
syst[lat.neighbors()] = -t

syst = kwant.wraparound.wraparound(syst)

leadL = kwant.Builder(kwant.TranslationalSymmetry((-1,0), lat.vec((0,W))))
leadL[(lat(0, j) for j in range(W))] = 4 * t
leadL[lat.neighbors()] = -t

syst.attach_lead(kwant.wraparound.wraparound(leadL, keep=0))

leadR = kwant.Builder(kwant.TranslationalSymmetry((1,0), lat.vec((0,W))))
leadR[(lat(0, j) for j in range(W))] = 4 * t
leadR[lat.neighbors()] = -t

syst.attach_lead(kwant.wraparound.wraparound(leadR, keep=0))

#kwant.plot(syst)
syst = syst.finalized()
#wf1 = kwant.solvers.default.wave_function(syst, energy=.2)(0)
plot_conductance(syst, energies=[.01 * i for i in range(200)])
