import kwant
from matplotlib import pyplot
# import scipy.sparse.linalg as sla
import numpy as np
import scipy.io as sio


def make_system(a=1, t=1.0, W=10, L=30):
    lat = kwant.lattice.square(a)
    syst = kwant.Builder()
    syst[(lat(x, y) for x in range(L) for y in range(W))] = 4 * t  # makes lattice
    syst[lat.neighbors()] = -t  # .neighbors adds hoppings

    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))  # Adding leads on left side
    lead[(lat(0, j) for j in range(W))] = 4 * t
    lead[lat.neighbors()] = -t

    lead2 = kwant.Builder(kwant.TranslationalSymmetry((a, 0)))  # Adding leads on left side
    lead2[(lat(0, j) for j in range(W))] = 4 * t
    lead2[lat.neighbors()] = -t
    # Attachs the leads onto the system. .reversed() adds the right leads
    syst.attach_lead(lead)
    syst.attach_lead(lead2)

    return syst


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


syst = make_system()

syst = syst.finalized()
# kwant.plot(syst)
# plot_conductance(syst, energies=[.01 * i for i in range(200)])
wf1 = kwant.solvers.default.wave_function(syst, energy=.7)(0)
#sio.savemat('wf/'+'wf6.mat', {'wf1':wf1})
kwant.plotter.map(syst, np.real(wf1[1]), cmap='jet')
