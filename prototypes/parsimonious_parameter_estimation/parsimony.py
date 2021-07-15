'''
Functions related to finding 'unit parsimonious perturbation vectors'.
'''

from __future__ import absolute_import, division, print_function

import numpy as np
from six.moves import range

def construct_parsimony_objects(
		perturbed_parameter_directions,
		restrained_parameter_directions
		):
	'''
	Computes and returns the objects needed for performing 'parsimonious'
	perturabtions.

	The math needed to arrive at these operations is long and requires some
	linear algebra background, although no one stage in the process is terribly
	sophisticated.  In short:

	There are two inputs; vectors describing directions in which we wish to
	perturb (encoded as rows in perturbed_parameter_directions), and vectors
	describing directions in which we wish to minimize changes (encoded as rows
	in the matrix restrained_parameter_directions).  In both cases, the
	elements of the matrices correspond to linear coefficients.

	For notational simplicity we'll consider one element (row) of
	perturbed_parameter_directions, which I will call 'd', and the full matrix
	restrained_parameter_directions, which I will call 'A'.  Then, a single
	unit parsimonious perturbation vector (a unit perturbation in the direction
	of 'd' while minimizing perturbations on 'A') is

	u = (I - N (A N)^p A) (d^T)^p

	I: identity matrix
	N: nullspace basis on d^T
	^p: Moore-Penrose pseudoinverse
	^T: transpose

	Two ways to obtain a nullspace basis are
	1) via portioning the matrices resulting from singular-value decomposition,
	2) by computing N = I - A^p A, where A is some matrix

	Herein I use the latter even though it results in a rank-deficient matrix.

	Parsimonious perturbation isn't a proper optimization strategy unless the
	perturbed_parameter_directions has full column rank.  An easy way to
	achieve this is to append an identity matrix to the bottom of the
	perturbed_parameter_directions matrix; however this results in several
	redundant 'u' vectors.  Consequently a post-processing step is to detect
	and discard these vectors; I likewise return the reduced
	perturbed_parameter_directions matrix so that the actual directions of
	perturbation are explictly encoded.
	'''

	(n_perturb, n_parameters) = perturbed_parameter_directions.shape

	n_parsimony = restrained_parameter_directions.shape[0]

	parsimonious_perturbation_directions = []

	for vector in perturbed_parameter_directions:
		# It's important that this be an explicit row vector
		vector_transpose = vector[None, :]
		vector_pseudoinverse = np.linalg.pinv(vector_transpose)

		# Nullspace projection matrix for the vector
		N = np.identity(n_parameters) - vector_pseudoinverse.dot(vector_transpose)

		# A second projection matrix (no known name, perhaps unnamed)
		projection_matrix = np.identity(n_parameters) - N.dot(
			np.linalg.pinv(restrained_parameter_directions.dot(N)).dot(
				restrained_parameter_directions
				)
			)

		parsimonious_perturbation_directions.append(
			projection_matrix.dot(vector_pseudoinverse)
			)

	parsimonious_perturbation_directions = np.concatenate(
		parsimonious_perturbation_directions,
		1
		).transpose()

	nonredundant = _nonredundant_directions(
		parsimonious_perturbation_directions,
		True
		)

	parsimonious_perturbation_directions = parsimonious_perturbation_directions[
		nonredundant, :
		]

	perturbed_parameter_directions = perturbed_parameter_directions[
		nonredundant, :
		]

	return parsimonious_perturbation_directions, perturbed_parameter_directions

_TOLERANCE = np.sqrt(np.finfo(np.float64).resolution)
def _nonredundant_directions(vectors, opposites_are_redundant = False):
	'''
	Returns a boolean vector that correspondes to a nonredundant set of
	vectors in the provided set of vectors.  Early vectors are favored over
	later vectors.
	'''

	vectors = np.array(vectors)

	n_vectors = vectors.shape[0]

	is_nonredundant = np.ones(n_vectors, bool)

	vectors = vectors / np.sqrt(np.sum(np.square(vectors), 1))[:, None]

	for i in range(n_vectors):
		vector_i = vectors[i, :]

		for j in range(i+1, n_vectors):
			if is_nonredundant[j]:
				vector_j = vectors[j, :]

				residual = np.dot(vector_i, vector_j)

				if opposites_are_redundant:
					residual = np.abs(residual)

				if np.abs(residual - 1.0) < _TOLERANCE:
					is_nonredundant[j] = False
					continue

	return is_nonredundant
