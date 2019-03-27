from sympy import Matrix, cos, sin
import numpy as np

def ode_func(par):

	return lambda t,y : [y[2], y[3],
 0.024*(30.0*y[1] - 0.006*y[3]**2*sin(y[1]) + 10.0*y[3] - 94.2477796076938)/(-3.6e-5*cos(y[1])**2 + 0.000744) - 0.0035316*sin(y[1])*cos(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744),
 0.006*(30.0*y[1] - 0.006*y[3]**2*sin(y[1]) + 10.0*y[3] - 94.2477796076938)*cos(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744) - 0.0182466*sin(y[1])/(-3.6e-5*cos(y[1])**2 + 0.000744)]