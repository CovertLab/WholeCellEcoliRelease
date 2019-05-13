"""
Utilities to compile functions, esp. from Sympy-constructed Matrix math.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from numba import njit
from sympy import Matrix
from typing import Callable


def build_function(arguments, expression, jit=True):
	# type: (str, str, bool) -> Callable
	"""Build a function from its arguments and source code expression, give it
	access to `Numpy as np`, and set up Numba to JIT-compile it on demand.

	Numba will optimize expressions like 1.0*y[2]**1.0 while compiling it
	to machine code.

	Args:
		arguments (str): comma-separated lambda argument names
		expression (str): expression to compile
		jit (bool): whether to JIT-compile the function; this option is just to
			work around Numba compiler bugs

	Returns:
		a Numba Dispatcher function(arguments)
	"""
	f = eval('lambda {}: {}'.format(arguments, expression), {'np': np}, {})

	if jit:
		# Too bad cache=True doesn't work with string source code.
		f_jit = njit(f, error_model='numpy')
		return f_jit
	return f


# TODO(jerry): Surely we can extract the argument array of "Matrix([...])" via
#  sympy calls more reliably than str(expr)[7:-1].
def _matrix_to_array(matrix):
	# type: (Matrix) -> str
	"""Convert a sympy Matrix expression to an 'np.array([...])' literal."""
	matrix_string = str(matrix)
	assert matrix_string.startswith('Matrix([')
	return 'np.array({})'.format(matrix_string[7:-1])


def derivatives(matrix, jit=True):
	# type: (Matrix, bool) -> Callable
	"""Build an optimized derivatives ODE function(y, t)."""
	return build_function('y, t',
		_matrix_to_array(matrix) + '.reshape(-1)', jit)

def derivatives_jacobian(jacobian_matrix, jit=True):
	# type: (Matrix, bool) -> Callable
	"""Build an optimized derivatives ODE Jacobian function(y, t)."""
	return build_function('y, t', _matrix_to_array(jacobian_matrix), jit)

def derivatives_with_rates(matrix, jit=True):
	# type: (Matrix, bool) -> Callable
	"""Build an optimized derivatives ODE function(y, t, kf, kr)."""
	return build_function('y, t, kf, kr',
		_matrix_to_array(matrix) + '.reshape(-1)', jit)

def derivatives_jacobian_with_rates(jacobian_matrix, jit=True):
	# type: (Matrix, bool) -> Callable
	"""Build an optimized derivatives ODE Jacobian function(y, t, kf, kr)."""
	return build_function('y, t, kf, kr',
		_matrix_to_array(jacobian_matrix), jit)
