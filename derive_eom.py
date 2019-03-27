#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 22:20:41 2019

@author: cldebuys
"""
import sys
sys.path.insert(0,"/home/cldebuys/Documents/Python/Tutorials") # this didnt work
import numpy as np
import sympy as sym
from multireplace import multireplace
from pythonFunction import pythonFunction
from sympy import Matrix, cos, sin
#import dill

# derives the eoms given the state pos, vel, and acc, and the Lagrangian L
def deriveEOM(q,dq,ddq,L):

    # x = np.concatenate((q,dq))
    x = Matrix([q,dq])

    # dx = np.concatenate((dq,ddq))
    dx = Matrix([dq,ddq])

    eom = (L.jacobian(dq)).jacobian(x)*dx - (L.jacobian(q)).reshape(len(q),1)
    
    # eom_full = np.concatenate((dq,eom))
    eom_full = Matrix([dq,eom])
    # returns the state vector, its derivative, and the equations of motion
    return(x,dx,eom)
   
def systemMatrices(ddq,eom):
    M = eom.jacobian(ddq)
    # H = C(q,dq)*dq + G(q)
    H = eom - M*ddq
    Minv = M.inv()
    return(M,H,Minv)
    
    
# =============================================================================
# # testing the above function    
# I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')
# 
# q0, q1, dq0, dq1, ddq0, ddq1, u = sym.symbols('q0, q1, dq0, dq1, ddq0, ddq1, u')
# 
# q =   sym.Matrix([  q0,   q1])
# dq =  sym.Matrix([ dq0,  dq1])
# ddq = sym.Matrix([ddq0, ddq1])
# 
# T = 0.5*(I0 + (m0+m1)*R**2)*dq[0]**2 - m1*l*R*dq[0]*dq[1]*sym.cos(q[1]) + 0.5*(I1 + m1*l**2)*dq[1]**2
# V = m1*g*l*(1 - sym.cos(q[1]))
# L = sym.Matrix([T - V])
# 
# x,dx,eom_full = deriveEOM(q,dq,ddq,L)
# 
# M,H,Minv = systemMatrices(ddq,eom_full)
# 
# # get eom for ode solver (eom_ode = ddq = Minv*H)
# eom_ode = Minv*(M*ddq - eom_full)
# 
# # sub in numbers for parameters
# par = (I0,I1,m0,m1,R,b,l)
# parval = (0.024,0.006,0.5,0.2,0.01,0.1,0.3)
# M_val = M.subs(zip(par, parval))
# 
# ##### GETTING EOM INTO ODE_FUNC #####
# # convert to string, edit values, and output to file
# M_str = str(M_val)
# 
# filename = "M_matrix.py" 
# funcname = "M_matrix"
# eqn = M_str
# varsToReplace = ["q0","q1","dq0","dq1"]
# importString = ["from sympy import Matrix, cos, sin"]
# ode_opt = 0
# pythonFunction(funcname,filename,eqn,varsToReplace,importString,[])
# 
# =============================================================================
