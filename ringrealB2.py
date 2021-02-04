import kwant                      # Recursive green function method 
import kwant.digest               # Random number generator built in kwant
import scipy.io as sio            # Module allowing to save in matlab format
import numpy.random               # Random number generator
import numpy as np                # Module with advanced math commands
from matplotlib import pyplot
#======================================================================
# Define the shape of the waveguide -----------------------------------
#======================================================================
def wv_shape(R):
    return ((R[0]**2 + R[1]**2) > (L-W/2)**2 and 
    (R[0]**2 + R[1]**2) < (L+W/2)**2)            # Const waveguide
#----------------------------------------------------------------------
a      = 1      # Lattice constant 
t      = 1      # Coupling between sites
E0     = 0      # On-site potential (-0.00029 - numerical and experimental)
d      = 0.16*t # Corresponds to experimental conditions: 0.5 -> ell=3.2um
L, W   = 1000, 120 
nrlz   = 5      # Number of realizations
E      = .35    # Energy N = ([W]'+1)/pi*acos((E-2)/2): 3.0-0.0003 for W=1000 3.0-0.0023 for W=2000 3.0-0.004 for W=3000 3.0-0.018 for W=150 3.0-0.009 for W=300 
phi0   = 0.00000
figsyn = 'n'    # 'y' to show figures and save them to files, 'n' - to save figures to files only
wfyn   = 'n'    # 'y' to compute wf, 'n' - do not
cvsb   = 'y'    # 'y' to compute conductance vs B flux
eigfun = 'n'    # 'y' to compute eigenfunctions
kxplot = 'n'    # 'y' to compute momenta for all of the modes vs B flux
GandN  = 'n'    # 'y' to compute conductance and modes for a system
tau_mat= 'y'    # 'y' to compute singular values for certain phis
if tau_mat == 'y':
    taus = np.empty((0, 0))
if cvsb == 'y':
    allcvsb = np.empty((0,0))
#------------------------------------------------------------------------------
def onsite(site, seed):
    return E0  + d  * (2.0*kwant.digest.uniform(repr(site),seed)-1.0)
def Aharonov_Bohm(site1, site2, phi):
    x1,y1=site1.pos
    x2,y2=site2.pos
    return np.exp(-1j * phi * (x1 - x2) * (y1 + y2) / np.power(L,2))  # phi is normalized for a period of 1
#------------------------------------------------------------------------------
lat  = kwant.lattice.square(a)
sys0 = kwant.Builder()
sys0[lat.shape(wv_shape, (L,0))] = onsite     # To make all sites the same
sys0[lat.neighbors()]           = Aharonov_Bohm     
#------------------------------------------------------------------------------
left_lead = kwant.Builder(kwant.TranslationalSymmetry((-1,0)))
left_lead[(lat(-L,y) for y in range(-int(W/2),int(W/2)))] = E0
left_lead[lat.neighbors()] = Aharonov_Bohm
sys0.attach_lead(left_lead)                          
#------------------------------------------------------------------------------
right_lead = kwant.Builder( kwant.TranslationalSymmetry( (1,0) ))
right_lead[(lat(L,y) for y in range(-int(W/2),int(W/2)))] = E0  
right_lead[lat.neighbors()] = Aharonov_Bohm
sys0.attach_lead(right_lead) 
#------------------------------------------------------------------------------
sys = sys0.finalized()
if(figsyn=='y'):
    kwant.plot(sys)   # Plots the shape of the waveguide
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
        p = dict(k_x=0,k_y=0,seed=str(rlz), phi = 0)
        s = kwant.smatrix(sys,E,params=p,check_hermiticity=True)
        wf0 = kwant.wave_function(sys,E,params=p,check_hermiticity=True)(0)
    #    sio.savemat('wf.mat', {'wf0':wf0})
    #--------- Plot conductance vs magnetic flux ---------------------------       
    if cvsb == 'y': 
        data = np.empty((0,0))
        phis = np.linspace(0,1,100)
        for phic in phis:
            p = dict(k_x=0,k_y=0,seed=str(rlz), phi=phic)
            smatrix = kwant.smatrix(sys,E,params=p)
            data = np.append(data, smatrix.transmission(1, 0))
        allcvsb = np.append(allcvsb, data)
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
    #------------------ Calculate Tau for many Phis ----------------------------
    if tau_mat == 'y':
        p       = dict(k_x=0,k_y=0,seed=str(rlz),phi=phi0)
        smatrix = kwant.smatrix(sys,E,params=p)
        t       = smatrix.submatrix(1,0)
        nL      = t.shape[1]
        steps = 100
        tau_phi = np.empty((0,0))
        phis = np.linspace(0,1,steps)
        print(phis)
        for phic in phis:
            p       = dict(k_x=0,k_y=0,seed=str(rlz),phi=phic)
            smatrix = kwant.smatrix(sys,E,params=p)
            t       = smatrix.submatrix(1,0)
            U,s,V  = np.linalg.svd(t,full_matrices=True,compute_uv=True) 
            tau    = np.power(np.abs(s[np.newaxis].T),2)
            print(tau.shape[0])
            if tau.shape[0] == 107:
                tau = np.append(tau, 0)
            tau_phi = np.append(tau_phi, tau)
        taus = np.append(taus, tau_phi)
    #--------------------Momentum plot---------------------------------------
    if kxplot == 'y':
        p = dict(k_x=0,k_y=0,seed=str(1),phi=.5)
        s = kwant.smatrix(sys,E,params=p,check_hermiticity=True)
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
if tau_mat == 'y':
    sio.savemat('tau_ring.mat', {'taus':taus, 'nrlz':nrlz, 'steps':steps,'nL':nL})
if cvsb == 'y':
    sio.savemat('cvsb_ring.mat', {'allcvsb':allcvsb})