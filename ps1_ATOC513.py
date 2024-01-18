"""

ps1_ATOC513.py

Author: Devon Gulley

Completion date: 2024/01/17

Sums together a number of 1-dimensional sinusoids to produce two different superposition patterns:
A Gaussian, and a wave packet both altered by a dispersion pattern omega in time.

Creates a small animation of each wave propagating in time.

Future intention is to render a 3d graph to show this evolution but this has yet to be implemented adequately.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

#Constants and problem setup
n_max = 100
epsilon = 0.025 #is very small
L = 100#due to the scaling of L and L_G, setting this will simply adjust the speed of the simulation
L_G = epsilon*L
A_0 = 2*epsilon*np.pi**0.5
x_ = np.linspace(start = -0.5*L, stop = 0.5*L, num = 1000)
g = 9.8
n0 = 1.5/epsilon#/L_G#3/epsilon #there must be a better way to do this, but 3/epsilon gave too many waves in the packet
display_height = 2 #from trial and error and is the vertical domain for which we can view the waves

#For repeated calls of A_n
def A(n):
    return A_0*np.exp(-(epsilon*n*np.pi)**2)

#For repeated calls of k_n
def k(n):
    return 2*n*np.pi/L


#eta(x,t) for the Gaussian
def gaussian(x,t):
    y = 0
    for n in range(1,n_max):
        omega = (g*k(n))**0.5
        y += A(n)*np.cos(k(n)*x - omega*t)
    return y

#eta(x,t) for the wave packet
def wavePacket(x,t):
    y = 0
    for n in range(1,n_max):
        B_n = A(n+n0)
        C_n = A(n-n0)
        omega = (g*k(n))**0.5

        y+= B_n*np.cos(k(n)*x - omega*t) + C_n*np.cos(k(n)*x - omega*t)

    return y



#Animated graphical visualisations of each wave
def wavePacketAnimation(t):
    ax.cla()
    y = wavePacket(x_,t)
    plt.plot(x_,y)
    view =display_height# 4*A_0;    
    ax.set_ylim(0-view,0+view)

def gaussianAnimation(t):
    ax1.cla()
    y = gaussian(x_,t)
    plt.plot(x_,y)
    view =display_height#4*A_0;    
    ax1.set_ylim(0-view,0+view)

#Displaying animated figures
fig = plt.figure()
ax = fig.add_subplot()
ani = animation.FuncAnimation(fig,wavePacketAnimation,frames= 600, interval=10)
plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot()
ani1 = animation.FuncAnimation(fig1,gaussianAnimation,frames= 600, interval=10)
plt.show()




"""
#There doesn't yet appear to be a readable way to plot this in 3D
#Here is some preliminary but unsatisfactory code


#graphing variables
xs = np.linspace(start = -0.5*L, stop = 0.5*L, num = 100)#number of x points to resolve to plot
ts = np.linspace(start = 0, stop = 30, num = 100)#length of time in seconds to plot

eta_packet = np.vectorize(wavePacket)
eta_gaussian = np.vectorize(gaussian)

fig2 = plt.figure()
ax2 = fig2.add_subplot(projection='3d')


import matplotlib.cm as cm
from matplotlib.colors import Normalize

cmap = cm.ocean
norm = Normalize(vmin=-1, vmax=1)


for t in ts:
    y = eta_gaussian(xs,t)
    ax2.scatter(xs,t ,y,color = cmap(norm(y)))
    ax2.set_xlabel('x')
    ax2.set_ylabel('t')
    ax2.set_zlabel('eta')
    ax2.set_title('Gaussian')
plt.show()
"""

