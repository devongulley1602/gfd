"""
bt.py

Author: Devon Gulley, Shayla Trembly

Completion date: 2024-04-18

The objective of this work is to illustrate some of the effects of Rossby waves under typical atmospheric conditions.

    - By implementing an equivalent barotropic PyQG model we can reasonably simulate 500mb mid-latitude atmospheric dynamics to first order
    - We use some initial conditions (namely Gaussian vortices) in the PV field to simulate a few examples of synoptic flow patterns and their
      relative progression across the spatial domain
    - Use simplified assumptions about the temporal and spatial scales
    - We impose on the output the generated velocity and height anomaly fields under shallow water potential equivalent barotropic assumptions
    - We constrain ourselves to the single layer horizontal case

This is to initialise the model, and is the framework for each initial condition case1.py to case5.py executed from the same directory.

"""

# Imports
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import pyqg
import multiprocessing
import os.path


# Simplifying assumptions about the time and space scales

# Scaling based on characteristic values from An Introduction to Dynamic Meteorology 2nd Edition - James Holton
# Chapter 4.5 A Scale Analysis of the Vorticity Equation
L = 1e6
H = 1e4
f0 = 1e-4
beta = 1e-11
U0 = 0.  # 1e-1#1e1#characteristic zonal wind


# Setting the characteristic values of the simulation and display
domain = 5  # how much wider the simulation is than the characteristic length
Lx = domain*L
Ly = domain*L
Ld = 0.1*L
Ro = Ld/L

# Creation of the spatial domain
nx = 512  # Grid point resolution
ny = nx
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)
# used for plotting vector fields (such as the wind) in a less cluttered way
spacing = 30*nx//256


T = 1e5  # characteristic time
tmax = T*10
dt = T*0.001

# Initialise the BTModel
threadsAvailable = multiprocessing.cpu_count()
m = pyqg.BTModel(L=Lx, W=Ly, nx=nx, ny=ny, H=H, beta=beta, rd=Ld,
                 U=U0, dt=dt, tmax=tmax, taveint=100*dt, ntd=threadsAvailable)


def gaussian_ellipse(x0=np.max(x)//6, y0=np.max(y)//2, sigma_x=Lx//64, sigma_y=Ly//64, intensity=1.0, velAttenuation=True):
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

    # Experimental
    attenuator_x = 1
    attenuator_y = 1
    if velAttenuation:
        # attenuates maximum velocities for different generated PV anomalies of this kind
        attenuator_x = (m.L//64)/sigma_x
        attenuator_y = (m.W//64)/sigma_y

    A = 21e-4  # A constant of proportionality
    d = np.sqrt((X-x0)**2 / (2*sigma_x)**2 + (Y-y0)**2 / (2*sigma_y) ** 2)
    qi = intensity * A * np.exp(-(d**2)) * \
        attenuator_x * attenuator_y
    return qi[np.newaxis, :, :]
