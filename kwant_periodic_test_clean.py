import kwant                      # Recursive green function method 
import kwant.digest               # Random number generator built in kwant
import scipy.io as sio            # Module allowing to save in matlab format
import numpy.random               # Random number generator
import numpy as np                # Module with advanced math commands

#======================================================================
# Define the shape of the waveguide -----------------------------------
#======================================================================
def wv_shape(pos):
    x, y = pos
    return (x>-1)&(x<L)&(y>-1)&(y<W1)   # Const waveguide
#----------------------------------------------------------------------
a      = 1     # Lattice constant 
t      = 1     # Coupling between sites
E0     = 0.0*complex(0,-0.0006) # On-site potential (-0.00029 - numerical and experimental)
E0L    = 0*t    # On-site potential in the lead
d      = 0.0*t # Corresponds to experimental conditions: 0.5 -> ell=3.2um
L      = 200   # Length of the waveguide 560->50um        - NEEDS TO BE EVEN
W1     = 100   # Width  of the waveguide 168->15um        - NEEDS TO BE EVEN
nrlz   = 1     # Number of realizations
E      = 3.749 # Energy N = ([W]'+1)/pi*acos((E-2)/2): 3.0-0.0003 for W=1000 3.0-0.0023 for W=2000 3.0-0.004 for W=3000 3.0-0.018 for W=150 3.0-0.009 for W=300 
figsyn = 'y'   # 'y' to show figures and save them to files, 'n' - to save figures to files only
wfyn   = 'y'   # 'y' to compute wf, 'n' - do not
dict_modes = {0: 0, 1: 1, 2: -1, 3: 2, 4: -2, 5: 3, 6: -3, 7: 4, 8: -4, 9: 5, 10: -5, 11: 6, 12: -6, 13: 7, 14: -7, 15: 8, 16: -8}
#------------------------------------------------------------------------------
def onsite(site,seed):
    return E0  + d  * (2.0*kwant.digest.uniform(repr(site),seed)-1.0)
# Define geometry -------------------------------------------------------------
lat  = kwant.lattice.square(a)
sys0 = kwant.Builder(kwant.TranslationalSymmetry( lat.vec((0,W1)) ))
#------------------------------------------------------------------------------
# Define onsite energies and couplings  ---------------------------------------
sys0[lat.shape(wv_shape,(0,0))] = onsite     # To make all sites the same
sys0[lat.neighbors()]           = t     
sys0 = kwant.wraparound.wraparound(sys0)
# -------------------------- Left lead ----------------------------------------
left_lead = kwant.Builder( kwant.TranslationalSymmetry((-1,0),lat.vec((0,W1))))
left_lead[(lat(0,y)   for y in range(W1))] = E0L
left_lead[lat.neighbors()] = t                                   
sys0.attach_lead(kwant.wraparound.wraparound(left_lead,keep=0)) 
# -------------------------- Right lead ---------------------------------------
right_lead = kwant.Builder( kwant.TranslationalSymmetry( (1,0), lat.vec((0,W1)) ))
right_lead[(lat(0,y) for y in range(W1))] = E0L    
right_lead[lat.neighbors()] = t                                  
sys0.attach_lead(kwant.wraparound.wraparound(right_lead,keep=0)) 
#------------------------------------------------------------------------------
sys = sys0.finalized()
if(figsyn=='y'):
    kwant.plot(sys)                                             # Plots the shape of the waveguide
#------------------------------------------------------------------------------
for rlz in range(nrlz):
    p = dict(k_x=0,k_y=0,seed=str(rlz))
    s = kwant.smatrix(sys,E,params=p,check_hermiticity=True)
    #--------------------------------------------------------------------------
    # On the first path save info about leads
    #--------------------------------------------------------------------------
    if (rlz==0):
        ypos=np.zeros((len(sys.leads[0].sites)))
        for j in range(len(sys.leads[0].sites)):
            ypos[j] = sys.leads[0].sites[j].pos[1]
        nL = s.submatrix(1,0).shape[1]
        #------ Define normalized profiles of the waveguide modes -------------
        khiL_n  = s.lead_info[0].wave_functions[:,:nL] # These are conjugated khi's
        khiL_n  = khiL_n/np.repeat( np.sqrt(np.sum( np.power(np.abs(khiL_n),2) ,axis=0))[np.newaxis,:], W1,axis=0)
        khiR_n  = np.zeros(khiL_n.shape)+complex(0,0)
        khiL0_n = np.zeros(khiL_n.shape)+complex(0,0)
        khiR0_n = np.zeros(khiL_n.shape)+complex(0,0)
        y=np.array(range(W1))
        for n in range(nL):
            khiL0_n[:,n] = np.exp( complex(0,1) * (2*np.pi*dict_modes[n]/W1) * y )
            khiL0_n[:,n] = khiL0_n[:,n] / np.power(np.dot( khiL0_n[:,n].conj().T , khiL0_n[:,n] ),0.5)
            khiR0_n[:,n] = khiL0_n[:,n] * np.exp( -complex(0,1) * (L-1) * np.arccos(E/2-np.cos(2*np.pi/W1*dict_modes[n])) ) # Mapping exp modes on Kwant modes, use: 
            khiR_n[ :,n] = khiL_n[ :,n] * np.exp( -complex(0,1) * (L-1) * s.lead_info[0].momenta[n]          ) # Mapping exp modes on Kwant modes, use: 
        #----- Propagation speed of modes in Kwant and exp bases --------------
        v                       = s.lead_info[0].velocities[nL:]
        v0                      = np.zeros((1,nL))
        v0[0,:int((nL-1)/2)+1 ] = s.lead_info[0].velocities[  nL  :    : 2]
        v0[0, int((nL-1)/2)+1:] = s.lead_info[0].velocities[2*nL-2:nL-1:-2]
        #------------- Define decomposition coefficients ----------------------
        projL = np.linalg.pinv( np.dot( khiL_n.T , khiL0_n.conj() ))
        projR = np.linalg.pinv( np.dot( khiR_n.T , khiR0_n.conj() ))
        sio.savemat('lead_info.mat', {'v':v,'v0':v0,'k':s.lead_info[0].momenta,'ypos':ypos,
                                      'psiL':s.lead_info[0].wave_functions,'psiL0':np.dot(projL,s.lead_info[0].wave_functions[:,:nL].T).T,
                                      'psiR':s.lead_info[1].wave_functions,'psiR0':np.dot(projR,s.lead_info[1].wave_functions[:,:nL].T).T,
                                      'khiL_n':khiL_n,'khiL0_n':khiL0_n,
                                      'khiR_n':khiR_n,'khiR0_n':khiR0_n, 'projL':projL,'projR':projR})
    #----------------- Compute wavefunctions ----------------------------------    
    if wfyn == 'y':
        wf0 = kwant.solvers.default.wave_function(sys,E,params=p,check_hermiticity=True)(0)
        wfe = np.dot( projL , wf0 ) # Wave functions in exp basis
        #----------------------------------------------------------------------
        t   = s.submatrix(1,0)  # This is incorrect
        #t  = np.dot( projR , np.dot( t , np.linalg.inv(projL) ) ) # This is incorrect because t is incorrect
        #----------------------------------------------------------------------
        tp  = np.dot(wf0[:,(L-1)*W1:],khiR_n.conj()).T * np.repeat(np.sqrt(v[:,np.newaxis]),nL,axis=1)
        te  = np.dot( projR , np.dot( tp, np.linalg.inv(projL) ) )
        #----------------------------------------------------------------------
        sio.savemat('wf.mat', {'wf0':wf0,'wfe':wfe,'t':t,'te':te})
        print( 'Calculation is done, g=',s.transmission(1,0))
    #--------------------------------------------------------------------------
#------------------------------------------------------------------------------
