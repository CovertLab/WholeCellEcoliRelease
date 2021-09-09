"""
Faster implementations of nonnegative least squares.
"""

import numpy as np
from scipy.optimize import nnls

def fast_nnls(A, b):
	"""
	Faster implementation of the nonnegative least squares algorithm, which
	returns a nonnegative vector x that minimizes ||Ax - b||_2. This function
	utilizes the property that both matrix A and vector b can be divided into
	matrices and vectors that each form a smaller nonnegative least squares
	problem, which can each be solved independently and the solutions later
	concatenated to yield the full vector x. Argument A is given as a sparse
	matrix.

	Args:
		A: scipy.sparse.csr.csr_matrix of size (M, N)
		b: numpy.ndarray of size (M, )
	Returns:
		x: numpy.ndarray of size (N, ), the solution to the NNLS problem.
		r: numpy.ndarray of size (M, ), the residual vector (Ax - b) of the NNLS
			problem.
	"""
	# Divide matrix A into smaller submatrices
	A_nonzero_row_indexes, A_nonzero_column_indexes = A.nonzero()

	visited_row_indexes = set()
	visited_column_indexes = set()
	submatrix_indexes = []

	def column_DFS(index, all_row_indexes, all_column_indexes):
		"""
		Recursive function to look for columns and rows in matrix A that should
		be grouped into the same NNLS problem.
		"""
		visited_column_indexes.add(index)
		all_column_indexes.append(index)

		for i in A_nonzero_row_indexes[A_nonzero_column_indexes == index]:
			if i not in visited_row_indexes:
				row_DFS(i, all_row_indexes, all_column_indexes)

	def row_DFS(index, all_row_indexes, all_column_indexes):
		"""
		Recursive function to look for columns and rows in matrix A that should
		be grouped into the same NNLS problem.
		"""
		visited_row_indexes.add(index)
		all_row_indexes.append(index)

		for i in A_nonzero_column_indexes[A_nonzero_row_indexes == index]:
			if i not in visited_column_indexes:
				column_DFS(i, all_row_indexes, all_column_indexes)

	# Loop through each column of matrix A
	for column_index in range(A_nonzero_column_indexes.max() + 1):
		# Search for columns and rows that can be grouped into a single NNLS
		# problem as the given column
		if column_index not in visited_column_indexes:
			submatrix_row_indexes = []
			submatrix_column_indexes = []
			column_DFS(column_index, submatrix_row_indexes, submatrix_column_indexes)

			submatrix_indexes.append((
				np.array(submatrix_row_indexes), np.array(submatrix_column_indexes)
				))

	# Initialize x
	x = np.zeros(A_nonzero_column_indexes.max() + 1)

	# Solve NNLS for each subproblem identified above
	for (row_indexes, column_indexes) in submatrix_indexes:
		if len(row_indexes) == 1 and len(column_indexes) == 1:
			x[column_indexes] = b[row_indexes]
		else:
			# Build a full submatrix A for each subproblem
			submatrix = np.zeros((len(row_indexes), len(column_indexes)))
			mask = np.isin(A_nonzero_row_indexes, row_indexes)
			for (i, j, v) in zip(
					A_nonzero_row_indexes[mask],
					A_nonzero_column_indexes[mask],
					A.data[mask]):
				submatrix[np.where(row_indexes == i)[0][0], np.where(column_indexes == j)[0][0]] = 1

			# Solve the subproblem
			x_subproblem, _ = nnls(submatrix, b[row_indexes])
			x[column_indexes] = x_subproblem

	assert np.all(x >= 0)

	# Calculate residuals vector
	r = A.dot(x) - b

	return x, r
