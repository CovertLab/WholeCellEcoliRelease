
from __future__ import division

import numpy as np

class RDPException(Exception): pass

def _squared_distance_to_line(first_point, last_point, points):
	"""
	Finds the distance between a line formed by two points and an
	additional set of points.

	This function uses vector projection and the Pythagorean theorem to
	find the orthogonal distance between a line and a point.  The
	orthogonal distance is the shortest possible distance from the point
	to a line.  We return the squared orthogonal distance to avoid a
	(costly) square root operation.

	Parameters
	----------

	first_point: A size m array, the first coordinate of the line.
	last_point: A size m array, the last coordinate of the line.

	points: A size n-by-m array, the coordinates of the points.

	Returns
	-------

	A size n array, the squared distances to the line.
	"""

	# TODO (John): explain/name intermediate calculations

	# Translate all points by the position of the first point
	last_point = last_point - first_point
	points = points - first_point

	# TODO (John): use a precision check rather than exact equality to zero
	if np.all(last_point == 0):
		# Degenerate case where the first and last point are identical,
		# collapsing the line to a point.  In this case, the distance
		# is just the distance from each point to the degenerate
		# endpoints.
		squared_distances = np.sum(points * points, 1)

	else:
		squared_distances = (
			np.sum(points * points, 1)
			- np.square(
				np.sum(points * last_point, 1)
				) / np.sum(last_point * last_point)
			)

	return squared_distances

def _fnz(vector):
	"""
	Returns the index of the first nonzero element in a vector.  No
	error checking.
	"""

	return np.where(vector)[0][0]

def rdp(points, threshold):
	"""
	Implementation of the Ramer-Douglas-Peucker (RDP) algorithm for
	reducing the number of points on a curve.

	The algorithm starts by forming a line between endpoints.  Then it
	calculates the Euclidean distance between that line and the
	intervening points.  If the distance falls below some threshold, the
	point is eliminated.  Otherwise the point becomes a new endpoint,
	and the procedure repeats until all points have been either set as
	endpoints or discarded.

	Parameters
	----------

	points: A size n-by-m array, where each row is a point on a curve,
	and each column is the location of that point in that dimension.
	Must be finite-valued and ordered.

	threshold: The maximum acceptable deviation from the line formed by
	two endpoints.

	Returns
	-------

	A size n bool array, true-valued wherever points should be retained.

	Notes
	-----

	General pattern of usage:

	>>> points = np.column_stack([x, y])

	>>> keep = rdp(points, threshold)

	>>> x_reduced = x[keep]
	>>> y_reduced = y[keep]

	For best results when plotting, the input data should be scaled
	according to the coordinates in figure space.  That is, any
	nonlinear scaling should be applied (i.e. log scaling), then the
	data should be normalized to the limits of the corresponding axis,
	and then finally, the data should be multiplied by the rendered size
	of the axis (in inches, for example).  Then the 'threshold'
	parameter can be interpreted in terms of physical units (again, in
	inches or similar).

	This algorithm could be implemented as a recursive function call;
	however, Python has a low default recursion limit.

	The 'threshold' parameter can be zero; this will only remove points
	that fall exactly on a line, up to limits of floating point
	precision.  Negative thresholds will be squared and therefore
	interpreted the same as corresponding positive thresholds.

	Sharp corners will generally be retained although there is no
	guarantee; the RDP algorithm does not make any assumptions about
	smoothness.

	If specific points along the curve must be retained, the coordinates
	should be split into groups at those points, passed individually to
	the RDP algorithm, and then the results stitched back together.

	The results of this function are deterministic and have no side
	effects; thus it is a true function.  Additionally, the result is
	unchanged if the order of the points is reversed except in
	degenerate cases.  The order of the dimensions does not matter.

	There are other Python implementations out there, but this is the
	fastest that I know of.  Cythonization would probably give minimal
	gains since the inner loop (distance calculations) is fully
	vectorized.
	"""

	points = np.asarray(points)

	if not np.all(np.isfinite(points)):
		raise RDPException('All points must be finite-valued.')

	(n_points, n_dims) = points.shape

	if (n_dims > n_points):
		print 'Number of dimensions appears to be greater than the number of elements; input may be transposed'

	# We work with squared distances to avoid calculating square roots, which is computationally expensive
	squared_threshold = np.square(threshold)

	keep = np.zeros(n_points, np.bool) # the positions to retain after filtering
	active = np.ones(n_points, np.bool) # points that have yet to be analyzed

	# Mark the first and last points as 1) kept and 2) analyzed
	keep[0] = True
	keep[-1] = True

	active[0] = False
	active[-1] = False

	first = 0

	# Repeat procedure until there are no more active points
	while np.any(active):
		# Find the next endpoints to test, starting from the left
		first = _fnz(~active[first:-1] & active[first+1:]) + first # The first endpoint is an inactive point just before an active point
		last = _fnz(~active[(first+1):]) + (first+1) # The last endpoint is the next inactive point after the first endpoint

		between = slice((first+1), last)

		square_distances = _squared_distance_to_line(
			points[first, :],
			points[last, :],
			points[between, :]
			)

		# TODO (John): catch and handle the degenerate cases where multiple points are at the maximum squared distance
		largest = np.argmax(square_distances)

		if square_distances[largest] > squared_threshold:
			# If the further point from the line formed by the endpoints
			# is more than a threshold, keep that point and use it as an
			# endpoint in subsequent iterations.
			keep[largest + (first+1)] = True
			active[largest + (first+1)] = False

		else:
			# Otherwise, discard all intervening points and move on.
			keep[between] = False
			active[between] = False

	return keep

def main():
	import matplotlib.pyplot as plt

	THRESHOLD = 1e-3 # in inches
	N = int(1e4)

	XLIM = (-0.5, 10.5)
	YLIM = (-0.9, +1.2)

	XSPAN = XLIM[1] - XLIM[0]
	YSPAN = YLIM[1] - YLIM[0]

	# These are just guesses at the size of the axes
	AXES_HEIGHT = 3 # in inches
	AXES_WIDTH = 4 # in inches

	test_x = np.linspace(0, 10, N)
	test_y = np.sin(test_x * np.pi) * np.exp(-test_x/2)

	points = np.column_stack([
		test_x/XSPAN*AXES_HEIGHT,
		test_y/YSPAN*AXES_WIDTH
		])

	keep = rdp(points, THRESHOLD)

	reduced_x = test_x[keep]
	reduced_y = test_y[keep]

	plt.subplot(1, 2, 1)

	plt.plot(test_x, test_y, '-', lw = 9, color = (0.8,)*3)

	plt.xlim(XLIM)
	plt.ylim(YLIM)

	plt.title('Original curve')

	plt.subplot(1, 2, 2)
	plt.plot(reduced_x, reduced_y, '-', lw = 9, color = (0.8,)*3)
	plt.plot(reduced_x, reduced_y, '.', ms = 9, color = (0.3,)*3)

	plt.title(
		'Retained {:0.2%} of points'.format(keep.sum()/keep.size)
		)

	plt.xlim(XLIM)
	plt.ylim(YLIM)

	plt.savefig('rdp_demo.png')

if __name__ == '__main__':
	main()
