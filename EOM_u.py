from sympy import Matrix, cos, sin
import numpy as np

def ode_func(t,y,par):
    
    Kp = 300.0
    Kd = 100.0
# =============================================================================
#     Kp = 30.0
#     Kd = 10.0
#     ref = np.pi
#     e = ref - q[1]
#     de = -q[3]
    #-(Kp*(np.polyval(par[1],t)-y[1])+Kd*(np.polyval(par[2],t)-y[3]))
    #-(Kp*(np.pi-y[1])-Kd*y[3])
# =============================================================================
    
    return [y[2], y[3],
 0.024*(-(Kp*(np.polyval(par[1],t)-y[1])+Kd*(np.polyval(par[2],t)-y[3])) - 0.006*y[3]**2*sin(y[1]))/(-3.6e-5*cos(y[1])**2 + 0.000744) - 0.0035316*sin(y[1])*cos(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744),
 0.006*(-(Kp*(np.polyval(par[1],t)-y[1])+Kd*(np.polyval(par[2],t)-y[3])) - 0.006*y[3]**2*sin(y[1]))*cos(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744) - 0.0182466*sin(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744)]
