"""

ps1_ATOC513.py

Author: Devon Gulley

Completion date: (ongoing)

Sums together a number of 1-dimensional sinusoids to produce two different superposition patterns:
A Gaussian, and a wave packet both altered by a dispersion pattern omega in time.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

#Constants and problem setup
n_max = 100
n0 = 60
epsilon = 0.025 #is very small
L = 1000
L_G = epsilon*L
A_0 = 2*epsilon*np.pi**0.5
x = np.linspace(start = -0.5*L, stop = 0.5*L, num = 1000)
g = 9.8

#For repeated calls of A_n
def A(n):
    return A_0*np.exp(-(epsilon*n*np.pi)**2)

#For repeated calls of k_n
def k(n):
    return 2*n*np.pi/L


#eta(x,y,t) for the Gaussian
def gaussian(t):
    y = 0
    for n in range(1,n_max):
        omega = (g*k(n))**0.5
        y += A(n)*np.cos(k(n)*x - omega*t)
    return y
    #plt.plot(x,y)
    #plt.show()

#eta(x,y,t) for the wave packet
def wavePacket(t):
    y = 0
    for n in range(1,n_max):
        B_n = A(n+n0)
        C_n = A(n-n0)
        omega = (g*k(n))**0.5

        y+= B_n*np.cos(k(n)*x - omega*t) + C_n*np.cos(k(n)*x - omega*t)

    return y

#y = gaussian(0)
#plt.plot(x,y)
#plt.show()

#y = wavePacket(0)
#plt.plot(x,y)
#plt.show()

def wavePacketAnimation(t):
    ax.cla()
    y = wavePacket(t)
    plt.plot(x,y)
    view =2# 4*A_0;    
    ax.set_ylim(0-view,0+view)

def gaussianAnimation(t):
    ax.cla()
    y = gaussian(t)
    plt.plot(x,y)
    view =2# 4*A_0;    
    ax.set_ylim(0-view,0+view)

fig = plt.figure()
ax = fig.add_subplot()
#ani = animation.FuncAnimation(fig,wavePacketAnimation,frames= 600, interval=10)
#plt.show()


ani = animation.FuncAnimation(fig,gaussianAnimation,frames= 600, interval=10)
plt.show()





