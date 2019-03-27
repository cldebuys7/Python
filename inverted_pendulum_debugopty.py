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
from opty.direct_collocation import Problem
from opty.utils import building_docs
import matplotlib.pyplot as plt
import matplotlib.animation as animation

target_angle = np.pi
duration = 1.3
num_nodes = 500
save_animation = False

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
I, m, g, d, t = sym.symbols('I, m, g, d, t')
theta, omega, T = sym.symbols('theta, omega, T', cls=sym.Function)

state_symbols = (theta(t), omega(t))
constant_symbols = (I, m, g, d)
specified_symbols = (T(t),)

eom = sym.Matrix([theta(t).diff() - omega(t),
                  I * omega(t).diff() + m * g * d * sym.sin(theta(t)) - T(t)])

# Specify the known system parameters.
par_map = OrderedDict()
par_map[I] = 1.0
par_map[m] = 1.0
par_map[g] = 9.81
par_map[d] = 1.0

# Specify the objective function and it's gradient.
def obj(free):
    """Minimize the sum of the squares of the control torque."""
    T = free[:1*num_nodes] - np.pi
    return interval_value * np.sum(T**2) # sum(T^2)
# =============================================================================
# def obj(free):
#     """Minimize the sum of the squares of the control torque."""
#     T = free[2 * num_nodes:]
#     return interval_value * np.sum(T**2) # sum(T^2)
# =============================================================================

def obj_grad(free):
    grad = np.zeros_like(free)
    grad[:1*num_nodes] = 2.0 * interval_value * (free[:1*num_nodes] - np.pi)
    return grad # 2T
# =============================================================================
# def obj_grad(free):
#     grad = np.zeros_like(free)
#     grad[2 * num_nodes:] = 2.0 * interval_value * free[2 * num_nodes:]
#     return grad # 2T
# =============================================================================

# Specify the symbolic instance constraints, i.e. initial and end
# conditions.
instance_constraints = (theta(0.0),
                        theta(duration) - target_angle,
                        omega(0.0),
                        omega(duration))

# Create an optimization problem.
prob = Problem(obj, obj_grad, eom, state_symbols, num_nodes, interval_value,
               known_parameter_map=par_map,
               instance_constraints=instance_constraints,
               bounds={T(t): (-10.0, 10.0)}) # (-2.0, 2.0)

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

# Make some plots
prob.plot_trajectories(solution)
prob.plot_constraint_violations(solution)
prob.plot_objective_value()

#plt.show()