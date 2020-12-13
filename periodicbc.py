import kwant
from matplotlib import pyplot
import numpy as np
import scipy.io as sio

def plot_conductance(syst, energies):
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy,params=dict(k_x = 0, k_y = 0))
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()

def make_system(a=1, t=1.0, L=100, W=48):
    lat = kwant.lattice.square(a)

    syst = kwant.Builder(kwant.TranslationalSymmetry(lat.vec((0,W))))
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

    syst = syst.finalized()
    return syst

syst = make_system()
#kwant.plot(syst)
wf0 = kwant.solvers.default.wave_function(syst, energy=.1, params=dict(k_x = 0, k_y = 0))(0)
#plot_conductance(syst, energies=[.04 * i for i in range(100)])
#kwant.plotter.map(syst, np.real(wf1[1]), cmap='jet')
sio.savemat('wf/scripts/'+'wfW36.mat', {'wf0':wf0})
