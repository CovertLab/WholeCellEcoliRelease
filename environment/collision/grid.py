from __future__ import absolute_import, division, print_function

import numpy as np

def raster(location, dx):
	''' translate a location from continuous values to an index into a grid of the given dx '''

	return map(int, np.floor(location / dx))

def within(bounds, point):
	''' determine whether the given point is within the provided bounding box '''

	return 0 <= point[0] < bounds[0] and 0 <= point[1] < bounds[1]


def normalize(vector):
	''' normalize the given vector based on its magnitude '''

	magnitude = np.linalg.norm(vector)
	normal = vector.copy()
	if magnitude == 0:
		normal *= 0
	else:
		normal /= magnitude

	return normal


class Shape(object):
	'''
	Base class for renderable shapes.

	This class provides a base class for defining shapes, and speifically the translation from
	some idealized expression of a shape into a series of indexes into a grid of a given dx.

	Clients of a Shape subclass will interact with it through the `indexes(dx)` method, which
	returns the rendered indexes if this dx has already been seen, or
	calls render if this dx is novel, memoizing the result.

	Subclassing Shape requires only defining the `__init__` method which accepts the information
	callers will be providing to express the shape (its minimal representation), and overriding the
	`render` method, which translates from the minimal represention to a series of indexes into a
	grid (the maximal representation). 
	'''

	def __init__(self):
		'''
		Setup the dictionary that will contain the various mappings from dx to indexes returned
		by `render(dx)`.
		'''
		self.renders = {}

	def center(self):
		'''
		Retrieve the center point of this shape.
		'''

		return np.array([0, 0])

	def indexes(self, dx):
		'''
		Return indexes if we have already calculated them for this dx, otherwise call render(dx)

		This is a memoization of `render(dx)`

		Args:
		    dx (float): resolution of the grid this shape will be applied to.

		Returns:
		    List[Tuple[int, int]]: a list of indexes into a grid of the given dx this shape occupies.
		'''

		if dx not in self.renders:
			self.renders[dx] = self.render(dx)
		return self.renders[dx]

	def render(self, dx):
		'''
		Translate from the minimal representation of this shape into a series of indexes into a
		grid of the given dx.

		Args:
		    dx (float): resolution of the grid this shape will be applied to.

		Returns:
		    List[Tuple[int, int]]: a list of indexes into a grid of the given dx this shape occupies.
		'''
		return []

class Line(Shape):
	'''
	Translate a begin and end point into a single index tuple for each row.
	'''

	def __init__(self, begin, end):
		'''
		Create a line with the given begin and end points.

		Args:
		    begin (Tuple[float, float]): A pair describing the begin point of the line.
		    end (Tuple[float, float]): A pair describing the end point of the line.
		'''
		super(Line, self).__init__()

		self.begin = np.array(begin)
		self.end = np.array(end)

	def render(self, dx):
		'''
		Render the line from begin and end points to indexes into a grid of the given dx.

		This function uses interpolation to find the points along the line. We first find the indexes
		for each row of the grid this line is present in, then for each index we find the ratio of
		how far along the line that index is and use this ratio to scale the begin and end
		points of the line. 

		Args:
		    dx (float): The resolution of the grid we are rendering on to.
		'''

		begin = raster(self.begin, dx)[1]
		end = raster(self.end, dx)[1] + 1
		interval = end - begin

		def interpolate(y):
			ratio = float(y - begin) / interval
			return self.begin * (1 - ratio) + self.end * ratio

		return [
			raster(interpolate(y), dx)
			for y in xrange(begin, end + 1)]

class Chain(Shape):
	'''
	Represent a sequence of lines as a series of endpoints.

	The Chain represents a sequence of lines as a series of points, with a line drawn between each
	point and the next point in the sequence. Once these are defined it uses a series of calls
	to Line to do the actual line rendering between the points.
	'''

	def __init__(self, points):
		'''
		Initialize a chain of lines with their endpoints.

		Args:
		    points (List[Tuple[float, float]]): A list of 2d points
		'''
		super(Chain, self).__init__()

		self.points = points

	def render(self, dx):
		'''
		Convert the list of points into a series of lines, then render those lines.

		Args:
		    dx (float): the incremental step size of the smallest unit within the grid being
		        rendered on to.
		'''

		cover = []
		if len(self.points) > 0:
			if len(self.points) > 1:
				first = Line(self.points[0], self.points[1])
				cover = first.render(dx)
				previous = self.points[1]
				for point in self.points[2:]:
					segment = Line(previous, point)
					previous = point
					squares = segment.render(dx)
					cover.extend(squares[1:])
			else:
				cover = [raster(self.points[0], dx)]

		return cover

class Rectangle(Shape):
	'''
	Translate a midpoint, a dimension, and an orientation into the 2d area it represents.

	The strategy for rendering a rectangle is to find the corners, render each side as a chain of
	lines (see `Chain` above), then using those lines as the bounds of a series of ranges, one
	for each row in the rectangle.
	'''

	def __init__(self, dimension, location, orientation):
		'''
		Construct a rectangle from its dimension, location (midpoint) and orientation.

		Args:
		    dimension (Tuple[float, float]): The width and height of this rectangle
		    location (Tuple[float, float]): The midpoint of this rectangle
		    orientation (float): The rotation angle of this rectangle
		'''

		super(Rectangle, self).__init__()

		self.dimension = np.array(dimension)
		self.location = np.array(location)
		self.orientation = orientation

	def center(self):
		return self.location

	def corners(self):
		'''
		Find the corners of this rectangle given its dimension, location and orientation

		Returns:
		    List[np.array[float]]: A four element list with one point for each corner
		'''

		cos = np.cos(self.orientation)
		sin = np.sin(self.orientation)
		diagonal = self.dimension * 0.5
		oneone = np.array([1, 1])

		return [
			self.location + (diagonal * quadrant).dot(oneone)
			for quadrant
			in np.array([
				[[-cos, -sin],
				 [-sin,  cos]],
				[[ cos, -sin],
				 [ sin,  cos]],
				[[-cos,  sin],
				 [-sin, -cos]],
				[[ cos,  sin],
				 [ sin, -cos]]])]

	def render(self, dx):
		'''
		Translate this rectangle's dimension, location, and orientation into a series of indexes
		into a grid with the given dx.

		Args:
		    dx (float): the smallest width represented in the grid.
		'''

		corners = self.corners()
		corners.sort(key=lambda x: x[1])
		begin = corners[0]
		a = corners[1]
		b = corners[2]
		end = corners[3]

		left = Chain([
			begin,
			a if a[0] < b[0] else b,
			end]).render(dx)
		right = Chain([
			begin,
			a if a[0] > b[0] else b,
			end]).render(dx)
		# perimeter = [left, right]

		indexes = np.concatenate([
			map(lambda x: [x, left[1]],
				range(left[0], right[0] + 1))
			for left, right
			in zip(left, right)])

		return set(map(tuple, indexes))

class Grid(object):
	'''
	A 2d grid of integer values representing the presence and overlap of shapes, useful for
	collision detection.

	The grid starts out empty (fill = -1) and accepts any number of `Shape` objects to impress
	into the grid. It achieves this by rendering the shapes into a series of indexes based on the
	dx of the grid, then incrementing the value at each index. Once this is complete any values
	in the grid greater than zero represent the overlap of two or more shapes. These values can then
	be read again by each shape to determine the total force on that shape (assuming each shape will
	be pushed away by overlap). The sum of all vectors from overlap to the midpoint of each shape
	represents the total force the shape is subjected to. 
	'''

	def __init__(self, bounds, dx):
		'''
		Create a grid out of the given bounds with the given dx.

		Args:
		    bounds (Tuple[float, float]): The outer bounds of this grid. It is assumed the origin of the
		        grid is (0, 0). TODO (Ryan): generalize this to a grid of any bounds.
		    dx (float): The resolution of the grid (the dimensions of the smallest area).
		'''

		self.bounds = np.array(bounds)
		self.dx = dx
		self.dimension = self.bounds / dx
		self.grid = np.full(map(int, self.dimension), -1)

	def reset(self):
		''' Reset the grid back to its empty state (fill = -1) '''

		self.grid.fill(-1)

	def impress(self, shape):
		'''
		Impress the area represented by this shape onto the grid. This is done by finding the
		indexes of the shape given this grid's dx, then incrementing the value at each index. Any
		indexes that lie outside the grid are ignored.

		Args:
		    shape (Shape): The shape we are going to impress.
		'''

		indexes = shape.indexes(self.dx)
		for index in indexes:
			if within(self.dimension, index):
				self.grid[index[0]][index[1]] += 1

	def overlap(self):
		''' Find the total overlap given by the previously impressed shapes '''

		total_overlap = np.sum(np.sum(self.grid[self.grid > 0]))
		return total_overlap

	def forces(self, shape):
		'''
		Given a shape, find the sum of all the vectors from any overlap to this shape's midpoint.

		To find the total force, we iterate through each index of the shape into this grid and find
		overlap, then find the vector from each overlapping index to this shape's `center()` and
		sum them. Finally, we divide by the square root of the total number of overlapping pixels
		to convert from 2d to 1d proportionality.

		Any index that is outside of the bounds of this grid is treated as overlap for the purposes
		of calculating this force.

		Args:
		    shape (Shape): The shape we are finding the forces for.
		'''

		location = np.array([int(i / self.dx) for i in shape.center()])
		total_force = np.array([0.0, 0.0])

		pixels = 0.
		for index in shape.indexes(self.dx):
			if within(self.dimension, index):
				# find collisions, get force vector for each.
				if self.grid[index] > 0:
					pixels += 1
					total_force += location - np.array(index)
			else:
				# check for boundaries
				if index[0] < 0:
					total_force[0] += 1
				if index[0] >= self.dimension[0]:
					total_force[0] -= 1
				if index[1] < 0:
					total_force[1] += 1
				if index[1] >= self.dimension[1]:
					total_force[1] -= 1
				pixels += 1

		# scaling = 15
		# total_force = normalize(total_force)

		scaling = 1.5 / np.sqrt(pixels) if pixels > 0 else 1
		total_force *= scaling * (self.dx ** 2)

		return total_force

	def collisions(self, shapes):
		'''
		Given a collection of shapes, find the total overlap and all forces on all shapes.

		Args:
		    shapes (list(Shape)): The shapes we will find the collisions for

		Returns:
		    overlap (int): Total overlap (value to be minimized)
		    forces (Dict[str, Tuple[float, float]]): A dictionary from shape id to total force vector for
		        that shape.
		'''

		self.reset()
		for shape_id, shape in shapes.iteritems():
			self.impress(shapes[shape_id])

		overlap = self.overlap()
		forces = {
			shape_id: self.forces(shape)
			for shape_id, shape in shapes.iteritems()}

		return overlap, forces
