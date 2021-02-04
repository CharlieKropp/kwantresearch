import kwant                      # Recursive green function method 
import kwant.digest               # Random number generator built in kwant
import scipy.io as sio            # Module allowing to save in matlab format
import numpy.random               # Random number generator
import numpy as np                # Module with advanced math commands
from matplotlib import pyplot
#======================================================================
# Define the shape of the waveguide -----------------------------------
#======================================================================
def wv_shape(pos):
    x, y = pos
    return (x>-1)&(x<L)&(y>-1)&(y<W1)   # Const waveguide
#----------------------------------------------------------------------
a      = 1      # Lattice constant 
t      = 1      # Coupling between sites
E0     = 0      # On-site potential (-0.00029 - numerical and experimental)
d      = 0.35*t  # Corresponds to experimental conditions: 0.5 -> ell=3.2um
L      = 600     # Length of the waveguide 560->50um        - NEEDS TO BE EVEN
W1     = 140     # Width  of the waveguide 168->15um        - NEEDS TO BE EVEN
nrlz   = 1      # Number of realizations
E      = .5   # Energy N = ([W]'+1)/pi*acos((E-2)/2): 3.0-0.0003 for W=1000 3.0-0.0023 for W=2000 3.0-0.004 for W=3000 3.0-0.018 for W=150 3.0-0.009 for W=300
phi0   = 0.00000
figsyn = 'n'    # 'y' to show figures and save them to files, 'n' - to save figures to files only
wfyn   = 'n'    # 'y' to compute wf, 'n' - do not
cvsb   = 'n'    # 'y' to compute conductance vs B flux
eigfun = 'y'    # 'y' to compute eigenfunctions
kxplot = 'n'    # 'y' to compute momenta for all of the modes vs B flux
GandN  = 'n'    # 'y' to compute conductance and modes for a system
#------------------------------------------------------------------------------
def onsite(site, seed):
    return E0  + d  * (2.0*kwant.digest.uniform(repr(site),seed)-1.0)
def Aharonov_Bohm(site1, site2, phi):
    return np.exp(-phi * complex(0,1) * (2 * np.pi))
#------------------------------------------------------------------------------
lat  = kwant.lattice.square(a)
sys0 = kwant.Builder(kwant.TranslationalSymmetry( lat.vec((0,W1)) ))
sys0[lat.shape(wv_shape,(0,0))] = onsite     # To make all sites the same
sys0[lat.neighbors()]           = t         
for hop in sys0.hoppings():
    if hop[0].tag[1]==11 and hop[1].tag[1]==10: 
        sys0[hop]=Aharonov_Bohm
sys0 = kwant.wraparound.wraparound(sys0)
#------------------------------------------------------------------------------
left_lead = kwant.Builder( kwant.TranslationalSymmetry((-1,0),lat.vec((0,W1))))
left_lead[(lat(0,y)   for y in range(W1))] = E0
left_lead[lat.neighbors()] = t
for hop in left_lead.hoppings():
    if hop[0].tag[1]==11 and hop[1].tag[1]==10: 
        left_lead[hop]=Aharonov_Bohm
sys0.attach_lead(kwant.wraparound.wraparound(left_lead,keep=0)) 
#------------------------------------------------------------------------------
right_lead = kwant.Builder( kwant.TranslationalSymmetry( (1,0), lat.vec((0,W1)) ))
right_lead[(lat(0,y) for y in range(W1))] = E0  
right_lead[lat.neighbors()] = t
for hop in right_lead.hoppings():
    if hop[0].tag[1]==11 and hop[1].tag[1]==10: 
        right_lead[hop]=Aharonov_Bohm
sys0.attach_lead(kwant.wraparound.wraparound(right_lead,keep=0)) 
#------------------------------------------------------------------------------
sys = sys0.finalized()
#if(figsyn=='y'):
#    kwant.plot(sys)   # Plots the shape of the waveguide
#------------------------------------------------------------------------------
for rlz in range(nrlz):
    #----------------- Compute Conductance and N-------------------------------
    if GandN == 'y':
        p = dict(k_x=0,k_y=0,seed=str(rlz), phi = 0)
        s = kwant.smatrix(sys,E,params=p,check_hermiticity=True)
        print( 'Number of modes, N =', s.submatrix(1,0).shape[0])
        print( 'Conductance, g =',s.transmission(1,0))  
    #----------------- Compute wavefunctions ----------------------------------    
    if wfyn == 'y':
        p = dict(k_x=0,k_y=0,seed=str(rlz), phi = .4)
        s = kwant.smatrix(sys,E,params=p,check_hermiticity=True)
        wf0 = kwant.wave_function(sys,E,params=p,check_hermiticity=True)(0)
        sio.savemat('wf.mat', {'wf0':wf0})
    #--------- Plot conductance vs magnetic flux ---------------------------       
    if cvsb == 'y': 
        data = []
        phis = np.linspace(0,2,50)
        for phic in phis:
            p = dict(k_x=0,k_y=0,seed=str(rlz), phi=phic)
            smatrix = kwant.smatrix(sys,E,params=p)
            data.append(smatrix.transmission(1, 0))
        pyplot.plot(phis, data)
        pyplot.xlabel('Phi')
        pyplot.ylabel('G')
        pyplot.title('Aharonov-Bohm Effect')
        pyplot.show()
    #-------------- Plot eigenfunctions -----------------------------------
    if eigfun == 'y':
        p       = dict(k_x=0,k_y=0,seed=str(rlz),phi=phi0)
        smatrix = kwant.smatrix(sys,E,params=p)
        t       = smatrix.submatrix(1,0)
        wf0     = kwant.solvers.default.wave_function(sys,E,params=p)(0)
        nL      = t.shape[1]
        U,s,V  = np.linalg.svd(t,full_matrices=True,compute_uv=True) 
        tau    = np.power(np.abs(s[np.newaxis].T),2)
        V      = np.conjugate(V.T)
        print(tau) 
        pyplot.ion()
        for n in range(nL):
            kwant.plotter.map(sys,np.real(np.dot( V[:,n].T, wf0[:] )), cmap='jet')
            print(np.sum(np.power(abs(np.dot(t,V[:,n])),2)))
        print('g=',smatrix.transmission(1,0))
        pyplot.pause(300)
    #--------------------Momentum plot---------------------------------------
    if kxplot == 'y':
        steps = 100
        nL0 = s.submatrix(1,0).shape[1]
        kx = np.empty((0, 0))
        phis = np.linspace(0, 1, steps)
        for i in range(steps):
            s = kwant.smatrix(sys,E,params=dict(k_x=0,k_y=0,seed=str(rlz), phi = phis[i]),check_hermiticity=True)
            nL = s.submatrix(1,0).shape[1]
            mode_k = s.lead_info[0].momenta[:nL]
            if nL != nL0:
                mode_k = np.insert(mode_k, 0, 0)
            for i in range(nL0):
                kx =np.append(kx, mode_k[i])
        kx = kx.reshape(steps,nL0)
        sio.savemat('kx.mat', {'kx':kx, 'phis': phis})
        pyplot.figure()
        for j in range(nL0):
            pyplot.plot(phis, kx[:, j])            
        pyplot.show()