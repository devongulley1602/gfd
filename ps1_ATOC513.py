import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#Constants and problem setup
epsilon = 0.025 #is very small
L = 1000
L_G = epsilon*L
A_0 = 2*epsilon*np.pi**0.5
y = 0
x = np.linspace(start = -0.5*L, stop = 0.5*L, num = 1000)

#Question 1
n_max = 100
for n in range(1,n_max):
    k = 2*n*np.pi/L
    A_n = A_0*np.exp(-(epsilon*n*np.pi)**2)
    y += A_n*np.cos(k*x)

#plt.plot(x,y)
#plt.show()


#Question 2
n_0 = 6
D = A_0*np.exp(-(epsilon*np.pi)**2)
g = 9.8

def eta(t):
    y = 0
    for n in range(1,n_max):
        B_n = D*np.exp(-(n-n_0)**2)
        C_n = D*np.exp(-(n+n_0)**2)
        k = 2*n*np.pi/L
        omega = (g*k)**0.5
        
        y += B_n*np.cos(k*x - omega*t)
    return y

y = eta(0)
plt.plot(x,y)
plt.show()




