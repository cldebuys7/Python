from sympy import Matrix, cos, sin

def M_matrix(y):

	return Matrix([[0.0240700000000000, -0.0006*cos(y[1])], [-0.0006*cos(y[1]), 0.0240000000000000]])
