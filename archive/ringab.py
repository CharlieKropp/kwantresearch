import numpy as np
from math import pi
import kwant
from matplotlib import pyplot
import scipy.io as sio            # Module allowing to save in matlab format

lat = kwant.lattice.square()
L,W=16,2
def myshape(R): return ( 
        (R[0]**2 + R[1]**2) > (L-W/2)**2 and                 
        (R[0]**2 + R[1]**2) < (L+W/2)**2)

H=kwant.Builder()

H[lat.shape(myshape,(L,0) )]=4

#H[lat.shape(myshape_ellipse,(int(L*1.14),0)) ]=4


H[lat.neighbors()]=1

def Aharonov_Bohm(site1,site2,phi): return np.exp(-2j*pi*phi)
    
for hop in H.hoppings():
    if hop[0].tag[0]==1 and hop[0].tag[1]>0 and hop[1].tag[0]==0: 
        H[hop]=Aharonov_Bohm

sym=kwant.TranslationalSymmetry(lat.vec((1,0)))
def lead_shape(R): return abs(R[1]) < W/2 and abs(R[0]) <3
Hlead =kwant.Builder(sym)
Hlead[lat.shape(lead_shape,(0,0) )]=4
Hlead[lat.neighbors()]=1
H.attach_lead(Hlead)
H.attach_lead(Hlead.reversed())

Hf=H.finalized()
wf0 = kwant.wave_function(Hf,3.3,params=dict(phi=.3),check_hermiticity=True)(0)
print(len(wf0))
kwant.plotter.map(Hf, np.real(wf0[0]), cmap='jet')

