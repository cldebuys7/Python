#!/usr/bin/env python2
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
#from opty.utils import state_derivatives
from opty.direct_collocation import Problem
#from opty.utils import building_docs
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation


target_position = 1.0
duration = 0.7
num_nodes = 20

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
# =============================================================================
# m, t = sym.symbols('m, t')
# x, xdot, u = sym.symbols('x, xdot, u', cls=sym.Function)
# =============================================================================

#
n = 1; # number of links
x = me.dynamicsymbols('x:{}'.format(n + 1))
xdot = me.dynamicsymbols('xdot:{}'.format(n + 1))
u = me.dynamicsymbols('u')

m = sym.symbols('m:{}'.format(n + 1))
l = sym.symbols('l:{}'.format(n))
I, g, t = sym.symbols('I, g, t')

state_symbols = (x, xdot)
constant_symbols = (I, m, l, g)
specified_symbols = (u,)
#

eom = sym.Matrix([x.diff() - xdot,
                  xdot.diff() - u])

# Specify the known system parameters.
par_map = OrderedDict()
par_map[I] = 1.0
par_map[m] = [1.0, 1.0]
par_map[g] = 9.81
par_map[l] = 1.0

# Specify the objective function and it's gradient.
def obj(free):
    """Minimize the sum of the squares of the control torque."""
    # Torque vector is last portion of solution vector
    u = free[2 * num_nodes:]
    return interval_value * np.sum(u**2) # sum(T^2)


def obj_grad(free):
    grad = np.zeros_like(free)
    grad[2 * num_nodes:] = 2.0 * interval_value * free[2 * num_nodes:]
    return grad # 2T

# Specify the symbolic instance constraints, i.e. initial and end
# conditions.
instance_constraints = (x(0.0),
                        x(duration) - target_position,
                        xdot(0.0),
                        xdot(duration))

# Create an optimization problem.  THIS CAUSES SPYDER ISSUES
prob = Problem(obj, obj_grad, eom, state_symbols, num_nodes, interval_value,
               instance_constraints=instance_constraints,
               bounds={u(t): (-10.0, 30.0)})

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

# Make some plots
prob.plot_trajectories(solution)
prob.plot_constraint_violations(solution)
prob.plot_objective_value()