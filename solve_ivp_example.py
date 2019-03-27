#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 12:59:11 2019

@author: cldebuys

This example is from:
 1  https://stackoverflow.com/questions/48428140/imitate-ode45-function-from-matlab-in-python
The second example is the python version of the matlab example from:
 2  https://www.mathworks.com/help/matlab/ref/ode45.html
"""

import numpy as np
import sympy as sym
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp  # for method 2

def vdp1(t, y):

    return np.array([y[1], (1 - y[0]**2)*y[1] - y[0]])


# method 1
t0, t1 = 0, 20                # start and end
t = np.linspace(t0, t1, 100)  # the points of evaluation of solution
y0 = [2, 0]                   # initial value
y = np.zeros((len(t), len(y0)))   # array for solution
y[0, :] = y0
r = integrate.ode(vdp1).set_integrator("dopri5")  # choice of method
r.set_initial_value(y0, t0)   # initial values
for i in range(1, t.size):
   y[i, :] = r.integrate(t[i]) # get one more value, add it to the array
   if not r.successful():
       raise RuntimeError("Could not integrate")
plt.plot(t, y)
plt.show()

# method 2
t0, t1 = 0, 20                # start and end
t = np.linspace(t0, t1, 100)  # the points of evaluation of solution
y0 = [2, 0]                   # initial value
sol = solve_ivp(vdp1, [t0, t1], y0, t_eval=t)
plt.plot(sol.t, sol.y.transpose())
plt.show()
plt.plot(sol.t, sol.y[0], sol.t, sol.y[1])
plt.show()