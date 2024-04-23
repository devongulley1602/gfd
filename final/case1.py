"""
case1.py

Small Gaussian experiencing mean flow

Author: Devon Gulley

Completion date: 2024-04-18
"""
from bt import *

print("Establishing initial conditions...")
m.U = 5.
q = np.zeros([1, ny, nx])
q += gaussian_ellipse(x0=m.L//3, y0=m.W//6, sigma_x=m.L // 10,
                      sigma_y=m.W//100, intensity=0.6)
m.set_q(q)
