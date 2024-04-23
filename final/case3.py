"""
case3.py

Large elliptical Gaussian wave train

Author: Devon Gulley

Completion date: 2024-04-18
"""
from bt import *

# Initial conditions
q = np.zeros([1, ny, nx])

numTrain = 4  # number of nodes in the wave train
sign = 1

for i in range(numTrain):
    q += sign*gaussian_ellipse(x0=i*m.L//numTrain + m.L//(2*numTrain),
                               y0=m.W//2, sigma_x=m.L//50, sigma_y=m.W//15, intensity=0.4)
    sign *= -1  # alternates PV sign to create positive-negative couplets
m.set_q(q)
