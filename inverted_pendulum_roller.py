# -*- coding: utf-8 -*-
"""This solves the simple pendulum swing up problem presented here:
http://hmc.csuohio.edu/resources/human-motion-seminar-jan-23-2014
A simple pendulum is controlled by a torque at its joint. The goal is to
swing the pendulum from its rest equilibrium to a target angle by minimizing
the energy used to do so.
"""

from collections import OrderedDict

import numpy as np
import sympy as sym
import sympy.physics.mechanics as me
from opty.direct_collocation import Problem
from opty.utils import building_docs
import matplotlib.pyplot as plt
import matplotlib.animation as animation

target_angle = np.pi
target_position = 0
duration = 0.5
num_nodes = 100
save_animation = False
n = 2 # number of degreees of freedom

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')

#q = me.dynamicsymbols('q:{}'.format(2*n))
q0, q1, q2, q3, u = sym.symbols('q0, q1, q2, q3, u', cls=sym.Function)
q = [q0(t), q1(t), q2(t), q3(t)]
state_symbols = q
#state_symbols = (q0(t), q1(t), q2(t), q3(t))
constant_symbols = (I0, I1, m0, m1, R, g, l) # missing b
specified_symbols = (u(t),)

# =============================================================================
# eom = sym.Matrix([q[0].diff() - q[2],
#                   (I0 + (m0 + m1)*R**2) * q[2].diff() + b*q[2] - m1*l*R * q[3].diff() * sym.cos(q[1]) + m1*l*R*q[3]**2*sym.sin(q[1]) - u(t),
#                   q[1].diff() - q[3],
#                   (I1 + m1*l**2) * q[3].diff() + m1 * g * l * sym.sin(q[1]) - m1*l*R * q[2].diff() * sym.cos(q[1]) + m1*l*R*q[2]*q[3]*sym.sin(q[1])])
# =============================================================================

# =============================================================================
eom = sym.Matrix([q[0].diff() - q[1],
                  (I0 + (m0 + m1)*R**2) * q[1].diff() + b*q[1] - m1*l*R * q[3].diff() * sym.cos(q[2]) + m1*l*R*q[3]**2*sym.sin(q[2]) - u(t),
                  q[2].diff() - q[3],
                  (I1 + m1*l**2) * q[3].diff() + m1 * g * l * sym.sin(q[2]) - m1*l*R * q[1].diff() * sym.cos(q[2]) + m1*l*R*q[1]*q[3]*sym.sin(q[2])])
# =============================================================================
print(eom)
temp0 = eom.diff(q[1].diff())
temp1 = eom.diff(q[3].diff())
M = sym.Matrix( [[temp0[1], temp1[1]], [temp0[3], temp1[3]]] )

# Specify the known system parameters.
par_map = OrderedDict()
par_map[I0] = 0.024
par_map[I1] = 0.006
par_map[m0] = 0.5
par_map[m1] = 0.2
# missing R
par_map[b] = 0.1
par_map[g] = 9.81
par_map[l] = 0.3

# Specify the objective function and it's gradient.


# =============================================================================
def obj(free):
    """Minimize the sum of the squares of the control torque."""
    u = free[4 * num_nodes:]
    cost = u**2
    return interval_value * np.sum(cost) # sum(T^2)
# =============================================================================

# =============================================================================
# def obj(free):
#     """Minimize the sum of the squares of the angle error."""
#     error = free[2 * num_nodes : 3 * num_nodes] - target_angle*np.ones_like(free[2 * num_nodes : 3 * num_nodes])
#     return interval_value * np.sum(error**2) # sum(T^2)
# =============================================================================

# =============================================================================
def obj_grad(free):
    grad = np.zeros_like(free)
    grad[4 * num_nodes:] = 2.0 * interval_value * free[4 * num_nodes:]
    return grad # 2T
# =============================================================================

# =============================================================================
# def obj_grad(free):
#     grad = np.zeros_like(free)
#     grad[4 * num_nodes:] = 2.0 * interval_value * free[2 * num_nodes : 3 * num_nodes]
#     return grad # 2T
# =============================================================================

# Specify the symbolic instance constraints, i.e. initial and end
# conditions.
    
# =============================================================================
# instance_constraints = (q0(0.0),
#                         q0(duration) - target_position,
#                         q1(0.0),
#                         q1(duration) - target_angle,
#                         q2(0.0),
#                         q2(duration),
#                         q3(0.0),
#                         q3(duration))
# =============================================================================

# =============================================================================
instance_constraints = (q0(0.0),
                        q0(duration) - target_position,
                        q1(0.0),
                        q1(duration),
                        q2(0.0),
                        q2(duration) - target_angle,
                        q3(0.0),
                        q3(duration))
# =============================================================================

# Create an optimization problem.
prob = Problem(obj, obj_grad, eom, state_symbols, num_nodes, interval_value,
               known_parameter_map=par_map,
               instance_constraints=instance_constraints,
               bounds={u(t): (-1000.0, 1000.0)})

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

# Make some plots
prob.plot_trajectories(solution)
#prob.plot_constraint_violations(solution)
prob.plot_objective_value()

# =============================================================================
# import sympy.utilities.lambdify as lambdify
# f = lambdify(((q[0], q[1], q[2], q[3], q[0].diff(), u(t)),), eom, modules=sym)
# =============================================================================

