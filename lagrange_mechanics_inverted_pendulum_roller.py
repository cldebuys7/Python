#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:18:03 2019

@author: cldebuys

This example is adapted from:
    https://docs.sympy.org/latest/modules/physics/mechanics/lagrange.html
"""

import sympy as sym
# from sympy.physics.mechanics import *
import sympy.physics.mechanics as me

I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')

q0, q2 = me.dynamicsymbols('q0 q2')
q0d, q2d = me.dynamicsymbols('q0 q2', 1) # the 1 is for the 1st derivative
T = 0.5*(I0 + (m0+m1)*R**2)*q0d**2 - m1*l*R*q0d*q2d*sym.cos(q2) + 0.5*(I1 + m1*l**2)*q2d**2
V = m1*g*l*(1 - sym.cos(q2))
L = T - V

LM = me.LagrangesMethod(L, [q0, q2])

me.mechanics_printing(pretty_print=False)
LM.form_lagranges_equations()

LM.mass_matrix
LM.forcing

print(LM.mass_matrix)
print(LM.forcing)

# with holonomic constraint
# =============================================================================
# LM = LagrangesMethod(L, [q1, q2], hol_coneqs=[q1 - q2])
# LM.form_lagranges_equations()
# # mass and force matrices augmented with constraints
# LM.mass_matrix_full
# LM.forcing_full
# print(LM.mass_matrix)
# print(LM.mass_matrix_full)
# =============================================================================

# If there are any non-conservative forces or moments acting on the system, 
# they must also be supplied as keyword arguments in a list of 2-tuples of the 
# form (Point, Vector) or (ReferenceFrame, Vector) where the Vector represents 
# the non-conservative forces and torques. Along with this 2-tuple, the 
# inertial frame must also be specified as a keyword argument. This is shown 
# below by modifying the example above
# =============================================================================
# N = ReferenceFrame('N')
# P = Point('P')
# P.set_vel(N, q1d * N.x)
# FL = [(P, 7 * N.x)]
# LM = LagrangesMethod(L, [q1, q2], forcelist=FL, frame=N)
# LM.form_lagranges_equations()
# =============================================================================

