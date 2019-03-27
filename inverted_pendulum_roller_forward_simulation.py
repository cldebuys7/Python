#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:03:46 2019

@author: cldebuys
"""
import numpy as np
import sympy as sym
from sympy import Matrix, sin, cos
from derive_eom import deriveEOM, systemMatrices
from pythonFunction import pythonFunction

I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')

q0, q1, dq0, dq1, ddq0, ddq1, u = sym.symbols('q0, q1, dq0, dq1, ddq0, ddq1, u')

q =   Matrix([  q0,   q1])
dq =  Matrix([ dq0,  dq1])
ddq = Matrix([ddq0, ddq1])

T = 0.5*(I0 + (m0+m1)*R**2)*dq[0]**2 - m1*l*R*dq[0]*dq[1]*cos(q[1]) + 0.5*(I1 + m1*l**2)*dq[1]**2
V = m1*g*l*(1 - cos(q[1]))
L = Matrix([T - V])

x,dx,eom = deriveEOM(q,dq,ddq,L)

M,H,Minv = systemMatrices(ddq,eom)

# get eom for ode solver (eom_ode = ddq = Minv*H)
Kp = 30.0
Kd = 10.0
ref = np.pi
e = ref - q[1]
de = -dq[1]
u = Kp*e + Kd*de
B = Matrix([1,0])
eom_ode = Minv*(M*ddq - eom - B*u)
eom_full = Matrix(np.concatenate((dq,eom_ode)))

# sub in numbers for parameters
par = (I0,I1,m0,m1,R,b,l,g)
parval = (0.024,0.006,0.5,0.2,0.1,0.1,0.3,9.81)
eom_val = eom_full.subs(zip(par,parval))
eom_val = np.squeeze(np.asarray(eom_val))
eom_str = np.array2string(eom_val, separator=', ')

filename = "EOM.py" 
funcname = "ode_func"
eqn = eom_str
varsToReplace = ["q0","q1","dq0","dq1"]
importString = ["from sympy import Matrix, cos, sin","import numpy as np"]
ode_opt = 1
pythonFunction(funcname,filename,eqn,varsToReplace,importString,ode_opt)

from EOM import ode_func
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

t0, t1 = 0, 3                 # start and end
t = np.linspace(t0, t1, 100)  # the points of evaluation of solution
y0 = [0, np.pi - np.pi/6 ,0,0]                   # initial value
sol = solve_ivp(ode_func, [t0, t1], y0, t_eval=t)
plt.plot(sol.t, sol.y.transpose())
plt.legend(['phi','theta','phidot','thetadot'])
plt.show()
plt.plot(sol.t, sol.y[1], sol.t, sol.y[3])
plt.legend(['theta','thetadot'])
plt.show()
