"""
Utilities to compile functions, esp. from Sympy-constructed Matrix math.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from numba import njit
from sympy import Matrix
from typing import Callable, Tuple


def build_functions(arguments, expression):
	# type: (str, str) -> Tuple[Callable, Callable]
	"""Build a function from its arguments and source code expression, give it
	access to `Numpy as np`, and set up Numba to JIT-compile it on demand.
	There will be overhead to compile the first time the jit version is called
	so two functions are returned and can be selected for optimal performance.

	Numba will optimize expressions like 1.0*y[2]**1.0 while compiling it
	to machine code.

	Args:
		arguments (str): comma-separated lambda argument names
		expression (str): expression to compile

	Returns:
		a lambda function(arguments)
		a Numba Dispatcher function(arguments)
	"""
	f = eval('lambda {}: {}'.format(arguments, expression), {'np': np}, {})

	# Too bad cache=True doesn't work with string source code.
	f_jit = njit(f, error_model='numpy')

	return f, f_jit


# TODO(jerry): Surely we can extract the argument array of "Matrix([...])" via
#  sympy calls more reliably than str(expr)[7:-1].
def _matrix_to_array(matrix):
	# type: (Matrix) -> str
	"""Convert a sympy Matrix expression to an 'np.array([...])' literal."""
	matrix_string = str(matrix)
	assert matrix_string.startswith('Matrix([')
	return 'np.array({})'.format(matrix_string[7:-1])


def derivatives(matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized derivatives ODE function(y, t)."""
	return build_functions('y, t',
		_matrix_to_array(matrix) + '.reshape(-1)')

def derivatives_jacobian(jacobian_matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized derivatives ODE Jacobian function(y, t)."""
	return build_functions('y, t', _matrix_to_array(jacobian_matrix))

def rates(matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized rates function(t, y, kf, kr)."""
	return build_functions('t, y, kf, kr',
		_matrix_to_array(matrix) + '.reshape(-1)')

def rates_jacobian(jacobian_matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized rates Jacobian function(t, y, kf, kr)."""
	return build_functions('t, y, kf, kr',
		_matrix_to_array(jacobian_matrix))
