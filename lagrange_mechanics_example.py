#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:18:03 2019

@author: cldebuys

This example is adapted from:
    https://docs.sympy.org/latest/modules/physics/mechanics/lagrange.html
"""

import sympy as sym
from sympy.physics.mechanics import *

I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')

q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', 1) # the 1 is for the 1st derivative
L = m0*q1d**2 + m1*q2d**2

LM = LagrangesMethod(L, [q1, q2])

mechanics_printing(pretty_print=False)
LM.form_lagranges_equations()

LM.mass_matrix
LM.forcing
print(LM.mass_matrix)

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

