'''
Linear algebra-related functions.
'''

from __future__ import absolute_import, division, print_function

import numpy as np
from six.moves import range

def approximate_gradient(f, x, dx):
	'''
	Computes the approximate gradient of a provided function f around a point x
	using center differences of sizes dx (single value, or vector).
	'''

	n = x.size

	df = np.empty(n, np.float64)

	if np.isscalar(dx):
		dx = dx * np.ones(n, np.float64)

	for i in range(n):
		xu = x.copy()
		xd = x.copy()

		xu[i] += dx[i]/2
		xd[i] -= dx[i]/2

		df[i] = (f(xu) - f(xd))/dx[i]

	return df

def approximate_gradient_adaptive(f, x, dx_initial, iteration_limit = 100, threshold = 1e-3, falloff = 0.1):
	'''
	Adaptively re-evaluates the gradient until it converges.
	'''

	old_gradient = None
	dx = dx_initial

	squared_threshold = threshold * threshold

	for i in range(iteration_limit):
		gradient = approximate_gradient(f, x, dx)

		if old_gradient is not None and np.sum(np.square(gradient - old_gradient)) < squared_threshold:
			break

		old_gradient = gradient
		dx *= falloff

	else:
		raise Exception('Gradient failed to converge')

	return gradient

