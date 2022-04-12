"""
Test the fast nonnegative least squares (NNLS) utils function

	cd wcEcoli
	pytest wholecell/tests/utils/test_fast_nnls.py
"""

from wholecell.utils.fast_nonnegative_least_squares import fast_nnls

import numpy as np
import numpy.testing as npt
from scipy import sparse
from scipy.optimize import nnls
import unittest
import time


def time_this(code_to_measure):
	"""
	Time the execution of code_to_measure() and return elapsed time in
	fractional seconds.
	"""
	elapsed_start = time.monotonic()
	code_to_measure()
	elapsed_end = time.monotonic()

	elapsed_time = elapsed_end - elapsed_start
	return elapsed_time


class Test_fast_nnls(unittest.TestCase):
	def setUp(self):
		self.default_array_size = 10
		np.random.seed(0)

	def test_return_value_dimensions(self):
		"""
		Test that return values have the correct dimensions.
		"""
		m = 5
		n = 3
		A = np.random.rand(m, n)
		b = np.random.rand(m)

		x, r = fast_nnls(A, b)
		assert x.shape == (n, )
		assert b.shape == (m, )

	def test_type_error(self):
		"""
		Test that arguments with wrong array types or dimensions raise
		TypeError exceptions.
		"""
		m = 5
		n = 3
		A = np.random.rand(m, n)
		sA = sparse.csr_matrix(A)
		A_wrongdim = np.random.rand(m)
		b = np.random.rand(m)
		sb = sparse.csr_matrix(b)
		b_wrongsize = np.random.rand(n)

		with self.assertRaisesRegex(TypeError, r'two-dimensional'):
			fast_nnls(A_wrongdim, b)
		with self.assertRaisesRegex(TypeError, r'one-dimensional'):
			fast_nnls(A, sb)
		with self.assertRaisesRegex(TypeError, r'one-dimensional'):
			fast_nnls(sA, sb)
		with self.assertRaisesRegex(TypeError, r'Dimensions of'):
			fast_nnls(sA, b_wrongsize)

	def test_identity_matrix(self):
		"""
		Test fast_nnls with an identity matrix A returns an array x that is
		equivalent to b.
		"""
		A = np.eye(self.default_array_size)
		b = np.random.rand(self.default_array_size)
		x, r = fast_nnls(A, b)

		npt.assert_array_equal(x, b)
		npt.assert_array_equal(r, np.zeros(self.default_array_size))

	def test_full_sparse_equivalence(self):
		"""
		Test if function returns same values for full and sparse matrix A's.
		"""
		A = np.random.rand(self.default_array_size, self.default_array_size)
		sA = sparse.csr_matrix(A)
		b = np.random.rand(self.default_array_size)

		x1, r1 = fast_nnls(A, b)
		x2, r2 = fast_nnls(sA, b)

		npt.assert_array_equal(x1, x2)
		npt.assert_array_equal(r1, r2)

	def test_reproducibility(self):
		"""
		Test reproducibility of fast_nnls outputs.
		"""
		A = np.random.rand(self.default_array_size, self.default_array_size)
		b = np.random.rand(self.default_array_size)

		x1, r1 = fast_nnls(A, b)
		x2, r2 = fast_nnls(A, b)

		npt.assert_array_equal(x1, x2)
		npt.assert_array_equal(r1, r2)

	def test_equilvalence_to_nnls(self):
		"""
		Test fast_nnls returns the same norm of the residual as scipy nnls for a
		random matrix A.
		"""
		A = np.random.rand(self.default_array_size, self.default_array_size)
		b = np.random.rand(self.default_array_size)

		_, rnorm_slow = nnls(A, b)
		_, r = fast_nnls(A, b)

		rnorm_fast = np.linalg.norm(r)

		self.assertAlmostEqual(rnorm_slow, rnorm_fast)

	def test_zero_column(self):
		"""
		Test fast_nnls returns a solution with a value of zero in the index
		corresponding to a column of zeros in matrix A.
		"""
		A = np.random.rand(self.default_array_size, self.default_array_size)
		b = np.random.rand(self.default_array_size)

		for i in range(A.shape[1]):
			A_copy = A.copy()
			A_copy[:, i] = 0
			x, _ = fast_nnls(A_copy, b)

			assert x[i] == 0

	def test_improved_performance(self):
		"""
		Test fast_nnls is faster than nnls for sparse arrays that can be
		decomposed into multiple nnls problems.
		"""
		# Decomposable array
		A = np.array([
			[1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
			])
		b = np.random.rand(10)

		time_slow = time_this(lambda: nnls(A, b))
		time_fast = time_this(lambda: fast_nnls(A, b))

		self.assertLess(time_slow, time_fast)
