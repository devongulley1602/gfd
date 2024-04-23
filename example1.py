"""
Example1.py

Author: Devon Gulley, Shayla Trembly 

Completion date: Ongoing

The objective of this work is to illustrate some of the effects of Rossby waves under typical atmospheric conditions.

    - By implementing an equivalent barotropic PyQG model we can reasonably simulate 500mb mid-latitude atmospheric dynamics to first order 
    - We use some initial conditions (namely Gaussian vortices) in the PV field to simulate a few examples of synoptic flow patterns and their 
      relative progression across the spatial domain
    - Use simplified assumptions about the temporal and spatial scales
    - We impose on the output the generated velocity and height anomaly fields under shallow water potential equivalent barotropic assumptions
    - We constrain ourselves to the single layer horizontal case

"""

#Imports
import matplotlib.pyplot as plt 
import numpy as np
from numpy import pi
import pyqg
import multiprocessing




#Simplifying assumptions about the time and space scales 


#Scaling based on characteristic values from An Introduction to Dynamic Meteorology 2nd Edition - James Holton 
#Chapter 4.5 A Scale Analysis of the Vorticity Equation
L = 1e6
H = 1e4
f0 = 1e-4
beta = 1e-11
U0 = 5.#1e-1#1e1 #characteristic zonal wind


#Setting the characteristic values of the simulation and display 
domain = 5 #how much wider the simulation is than the characteristic length
Lx = domain*L
Ly = domain*L
Ld = 0.1*L
#Creation of the spatial domain
nx = 256 #Grid point resolution
ny = nx
x = np.linspace(0,Lx,nx)
y = np.linspace(0,Ly,ny)
X,Y = np.meshgrid(x,y)
spacing = 10*nx//256 #used for plotting vector fields (such as the wind) in a less cluttered way


T = 1e5 #characteristic time
tmax = T*10
dt = T*0.001

#Initialise the BTModel
threadsAvailable = multiprocessing.cpu_count()
m = pyqg.BTModel(L = Lx, W = Ly,nx = nx, ny = ny,H = H, beta = beta , rd = Ld,U = U0, dt = dt,tmax = tmax , taveint = 100*dt,ntd = threadsAvailable)


def gaussian_storm(x0=np.max(x)//64, y0 = np.max(y)//2, sigma_x = Lx//64,sigma_y = Ly//64,intensity = 1.0,velAttenuation = True,spreadFactor = 1.):
    """

    Creates a Gaussian-shaped PV anomaly scaled to generate roughly the same vorticity.
    Enable velAttenuation to somewhat constrain maximum winds despite different spatial coverage.

    Parameters
    ----------
    x0,y0 : float
        The position of where the anomaly is centred

    sigma_x,sigma_y : float
        Proportional to the overall position spread

    intensity : float
        An arbitrary number from 0 to 1 to produce the maximum vorticity generated suitable for this example problem

    velAttenuation : bool
        Incomplete, true when to constrain maximum wind speeds for vortices with different spread

    """

    #Experimental
    attenuator_x = 1
    attenuator_y = 1
    if velAttenuation:
        attenuator_x = (m.L//64)/sigma_x # attenuates maximum velocities for different generated PV anomalies of this kind 
        attenuator_y = (m.W//64)/sigma_y

    A = 21e-4 #A constant of proportionality
    d = np.sqrt((X-x0)**2 /(2*sigma_x)**2 + (Y-y0)**2 / (2*sigma_y) **2)
    qi = intensity* A * np.exp(-(d**2)/spreadFactor) * attenuator_x * attenuator_y
    return qi[np.newaxis,:,:]


def create_gaussian(m=m,x0 = m.L//2, y0 = m.W//2,R = m.L//24):
    # # Gaussian IC
    fk = m.wv != 0
    ckappa = np.zeros_like(m.wv2)
    ckappa[fk] = np.sqrt( m.wv2[fk]*(1. + (m.wv2[fk])**2) )**-1
    nhx,nhy = m.wv2.shape
    Pi = -np.exp(-((m.x-x0)**2 + (m.y-y0)**2)/R**2)
    Pi = Pi - Pi.mean()
    Pi_hat = m.fft( Pi[np.newaxis,:,:] )
    KEaux = m.spec_var(m.wv*Pi_hat )
    pih = ( Pi_hat/np.sqrt(KEaux) )
    qih = -m.wv2*pih
    qi = m.ifft(qih)
    return qi

###############################################################################################################
#
#       
#       Specifying the initial conditions...
#
#
###############################################################################################################

q = np.zeros([1,ny,nx])
#q= gaussian_storm(x0=m.L//4,y0=m.W//4,sigma_x=m.L//100,sigma_y=m.W//100,intensity = 1,velAttenuation=False)
q= gaussian_storm(x0=m.L//3,y0=m.W//6,sigma_x=m.L//10,sigma_y=m.W//100,intensity = 0.2,velAttenuation=True)
q+= gaussian_storm(x0=m.L//6,y0=m.W//4,sigma_x=m.L//10,sigma_y=m.W//100,intensity = 0.6,velAttenuation=True)

#q = 3*create_gaussian(x0=m.L//5,y0 =m.W//4)
#q+= gaussian_storm(x0 = m.L//4,y0 = m.W*3/4,intensity = 0.1)
# np.max(q[0])
# q = create_gaussian(m)
# q += create_gaussian(m, x0 = m.L/4)
# q += create_gaussian(m, x0 = m.L*3/4)

m.set_q(q)

###############################################################################################################
#
#       
#
#
###############################################################################################################
input("...")
#Refresh and run the model
#refreshModel()
picNum = 0
for snapshot in m.run_with_snapshots(tsnapstart=0, tsnapint=50*m.dt):
    plt.clf()
    plt.title(f'Days ~{np.round(m.t/86400,2)}')
    picNum += 1
    U = (m.u[0] + m.U)*1.94
    V = (m.v[0])*1.94
    plt.barbs(X[::spacing, ::spacing], Y[::spacing, ::spacing], U[::spacing, ::spacing], V[::spacing, ::spacing], length=4, barbcolor='k', flagcolor='r', linewidth=0.4, sizes=dict(emptybarb=0.00))

    plt.imshow(m.q.squeeze(),origin='lower',extent = [0,m.L,0,m.W])
    
    plt.pause(1/60)


    plt.draw()
    plt.savefig(f'snapshots/interesting/{picNum}.jpg')
    plt.ioff()
