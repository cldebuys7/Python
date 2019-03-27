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
target_position = 0
duration = 0.5
num_nodes = 100
save_animation = False

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
I, M, m, b, g, l, t = sym.symbols('I, M, m, b, g, l, t')
x, xdot, theta, omega, u = sym.symbols('x, xdot, theta, omega, u', cls=sym.Function)

test = [x(t), xdot(t), theta(t), omega(t)]
state_symbols = (test)
constant_symbols = (I, M, m, g, l)
specified_symbols = (u(t),)

eom = sym.Matrix([x(t).diff() - xdot(t),
                  (M+m) * xdot(t).diff() + b*xdot(t) + m*l * omega(t).diff() * sym.cos(theta(t)) - m*l*omega(t)**2*sym.sin(theta(t)) - u(t),
                  theta(t).diff() - omega(t),
                  (I + m*l**2) * omega(t).diff() + m * g * l * sym.sin(theta(t)) + m*l * xdot(t).diff() * sym.cos(theta(t))])

# Specify the known system parameters.
par_map = OrderedDict()
par_map[I] = 0.006
par_map[M] = 0.5
par_map[m] = 0.2
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
instance_constraints = (x(0.0),
                        x(duration) - target_position,
                        xdot(0.0),
                        xdot(duration),
                        theta(0.0),
                        theta(duration) - target_angle,
                        omega(0.0),
                        omega(duration))

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
# # Display animation
# if not building_docs():
#     time = np.linspace(0.0, duration, num=num_nodes)
#     angle = solution[:num_nodes]
# 
#     fig = plt.figure()
#     ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-2, 2),
#                          ylim=(-2, 2))
#     ax.grid()
# 
#     line, = ax.plot([], [], 'o-', lw=2)
#     time_template = 'time = {:0.1f}s'
#     time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
# 
#     def init():
#         line.set_data([], [])
#         time_text.set_text('')
#         return line, time_text
# 
#     def animate(i):
#         x = [0, par_map[l] * np.sin(angle[i])]
#         y = [0, -par_map[l] * np.cos(angle[i])]
# 
#         line.set_data(x, y)
#         time_text.set_text(time_template.format(i * interval_value))
#         return line, time_text
# 
#     ani = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
#                                   interval=25, blit=True, init_func=init)
# 
#     if save_animation:
#         ani.save('pendulum_swing_up.mp4', writer='ffmpeg',
#                  fps=1 / interval_value)
# 
# plt.show()
# =============================================================================
