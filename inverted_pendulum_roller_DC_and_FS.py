#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
""" This code is for a wheeled inverted pendulum. The equations of motion are
    derived symbolically. A trajectory is generated
    for the controller using opty, and a forward simulation is conducted using
    solve_ivp.
"""

from collections import OrderedDict

import numpy as np
import sympy as sym
from sympy import Matrix, sin, cos
from derive_eom import deriveEOM, systemMatrices
import sympy.physics.mechanics as me
from opty.direct_collocation import Problem
from pythonFunction import pythonFunction

initial_angle = np.pi - np.pi/6
target_angle = np.pi
target_position = 0
duration = 1.0
num_nodes = 100
save_animation = False
n = 2 # number of degreees of freedom

interval_value = duration / (num_nodes - 1)

# Define symbols for system parameters
I0, I1, m0, m1, R, b, g, l, t = sym.symbols('I0, I1, m0, m1, R, b, g, l, t')
# Define the state vector
q = me.dynamicsymbols('q:{}'.format(2*n))
q0, q1, q2, q3, u = sym.symbols('q0, q1, q2, q3, u', cls=sym.Function)
u = sym.symbols('u', cls=sym.Function)
#q = [q0(t), q1(t), q2(t), q3(t)]
state_symbols = q
specified_symbols = (u(t),)

Pos =   Matrix([ q[0],  q[1]])
dPos =  Matrix([ q[2],  q[3]])
ddPos = Matrix([ q[2].diff(), q[3].diff()])
T = 0.5*(I0 + (m0+m1)*R**2)*q[2]**2 - m1*l*R*q[2]*q[3]*cos(q[1]) + 0.5*(I1 + m1*l**2)*q[3]**2
V = m1*g*l*(1 - cos(q[1]))
L = Matrix([T - V])
#eqn = (L.jacobian(dq)).jacobian(x)*dx - (L.jacobian(q)).reshape(len(q),1)
x,dx,eqn = deriveEOM(Pos,dPos,ddPos,L)
extra = Matrix([b*q[2]-u(t),0])
eom2 = eqn + extra
eom4 = Matrix([-dPos + Pos.diff(),eom2])

""" DIRECT COLLOCATION BEGIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> """
# Specify the known system parameters.
par_map = OrderedDict()
par_map[I0] = 0.024
par_map[I1] = 0.006
par_map[m0] = 0.5
par_map[m1] = 0.2
par_map[R] = 0.1
par_map[b] = 0.1
par_map[g] = 9.81
par_map[l] = 0.3

# Specify the objective function and it's gradient.

# =============================================================================
def obj(free):
    """Minimize the sum of the squares of the control torque."""
    #R = 10
    #e = (np.pi*np.ones_like(free[num_nodes:2*num_nodes]) - free[num_nodes:2*num_nodes])*R
    #e = (target_angle - free[num_nodes:2*num_nodes])*R
    u = free[4 * num_nodes:]
    #cost = e**2 + u**2
    cost = u**2
    return interval_value * np.sum(cost) # sum(T^2)
# =============================================================================

# =============================================================================
def obj_grad(free):
    #R = 10
    grad = np.zeros_like(free)
    grad[4 * num_nodes:] = 2.0 * interval_value * free[4 * num_nodes:]
    #grad[num_nodes:2*num_nodes] = 2.0 * interval_value * R*(np.pi*np.ones_like(free[num_nodes:2*num_nodes]) - free[num_nodes:2*num_nodes])
    return grad # 2T
# =============================================================================

# Specify the symbolic instance constraints, i.e. initial and end
# conditions.
# =============================================================================
instance_constraints = (q0(0.0),
                        q0(duration) - target_position,
                        q1(0.0) - initial_angle,
                        q1(duration) - target_angle,
                        q2(0.0),
                        q2(duration),
                        q3(0.0),
                        q3(duration))
# =============================================================================

# Create an optimization problem.
prob = Problem(obj, obj_grad, eom4, state_symbols, num_nodes, interval_value,
               known_parameter_map=par_map,
               instance_constraints=instance_constraints,
               bounds={u(t): (-1000.0, 1000.0)})

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

time = np.linspace(0, duration, num_nodes)
opt_control = solution[num_nodes*4:]
control_poly = np.polyfit(time,opt_control,16)
opt_angle = solution[num_nodes:num_nodes*2]
opt_vel = solution[num_nodes*3:num_nodes*4]
angle_poly = np.polyfit(time,opt_angle,16)
vel_poly = np.polyfit(time,opt_vel,16)
# Make some plots
prob.plot_trajectories(solution)
#prob.plot_constraint_violations(solution)
prob.plot_objective_value()
""" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DIRECT COLLOCATION END """

""" FORWARD SIMULATION BEGIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> """

# get eom for ode solver (eom_ode = ddq = Minv*H)
dq = dPos
ddq = ddPos
M,H,Minv = systemMatrices(ddq,eqn)

Kp = 30.0
Kd = 10.0
ref = np.pi
e = ref - q[1]
de = -q[3]
#u = Kp*e + Kd*de
u = sym.symbols('u')
B = Matrix([1,0])
eom_ode = Minv*(M*ddq - eqn - B*u)
eom_full = Matrix(np.concatenate((dq,eom_ode)))

# sub in numbers for parameters
par = (I0,I1,m0,m1,R,b,l,g)
parval = (0.024,0.006,0.5,0.2,0.1,0.1,0.3,9.81)
eom_val = eom_full.subs(zip(par,parval))
eom_val = np.squeeze(np.asarray(eom_val))
eom_str = np.array2string(eom_val, separator=', ')

# =============================================================================
# # generate function for ode to call (may need to edit after generation)
# filename = "EOM.py" 
# funcname = "ode_func"
# eqn2 = eom_str
# varsToReplace = ["q0(t)","q1(t)","q2(t)","q3(t)"]
# importString = ["from sympy import Matrix, cos, sin","import numpy as np"]
# ode_opt = 1
# pythonFunction(funcname,filename,eqn2,varsToReplace,importString,ode_opt)
# =============================================================================

from EOM_u import ode_func
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

t0, t1 = 0, 1.0                 # start and end
t_eval = np.linspace(t0, t1, 100)  # the points of evaluation of solution
y0 = [0, np.pi - np.pi/6 ,0,0]                   # initial value
y = np.zeros((len(y0),len(t_eval)))
# sol=solve_ivp(fun=lambda t, y: fun(t, y, *args), ...)
REF = [control_poly,angle_poly,vel_poly]
ref2 = 0

def wrapper(t, y):
    ode_func(t,y)
#sol = solve_ivp(lambda t,y: ode_func(t,y,ref2), [t0, t1], y0, t_eval=t_eval)
sol = solve_ivp(lambda t,y: ode_func(t,y,REF), t_span=[t0, t1], y0=y0, t_eval=t_eval)
# forward simulation trajectories
plt.figure(3)
plt.plot(sol.t, sol.y.transpose())
plt.legend(['phi','theta','phidot','thetadot'])
plt.show()

# forward simulation theta and thetadot
plt.figure(4)
plt.plot(sol.t, sol.y[1],'C1', sol.t, sol.y[3],'C3')
plt.legend(['theta','thetadot'])
plt.show()

# check polynomial fit of control input
plt.figure(5)
control_eval = np.polyval(control_poly,time)
plt.plot(time, opt_control, time, control_eval)
plt.legend(['optimal control','polynomial fit'])
plt.show()

# prepare optimal solution for plotting
opt_y0 = solution[:num_nodes*1]
opt_y1 = solution[num_nodes*1:num_nodes*2]
opt_y2 = solution[num_nodes*2:num_nodes*3]
opt_y3 = solution[num_nodes*3:num_nodes*4]
opt_u = solution[num_nodes*4:]
opt_y = [opt_y0,opt_y1,opt_y2,opt_y3]
opt_y = np.array(opt_y)

# compare theta & thetadot from forward sim to optimal trajectories
plt.figure(6)
plt.plot(sol.t, sol.y[1],'C1', time, opt_y1,'k--',sol.t, sol.y[3],'C3', time, opt_y3,'k--')
plt.legend(['forward sim theta','optimal traj','forward sim thetadot','optimal traj'])
plt.show()

# compare phi & phidot from forward sim to optimal trajectories
plt.figure(7)
plt.plot(sol.t, sol.y[0],'C0', time, opt_y0,'k--',sol.t, sol.y[2],'C2', time, opt_y2,'k--')
plt.legend(['forward sim phi','optimal traj','forward sim phidot','optimal traj'])
plt.show()

# actual contorl input
Kp = 300.0
Kd = 100.0
cnt = 0
actual_u = np.zeros_like(sol.t)
for i in sol.t:
    actual_u[cnt] = -(Kp*(np.polyval(REF[1],i)-sol.y[1][cnt])+Kd*(np.polyval(REF[2],i)-sol.y[3][cnt]))
    cnt = cnt + 1

plt.figure(8)
plt.plot(sol.t, actual_u,'b', time, opt_u,'k--')
plt.legend(['forward sim control','optimal control'])
plt.show()

""" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FORWARD SIMULATION END """

# dont use random initial guess
# debug why position error isnt minimized at all
# 