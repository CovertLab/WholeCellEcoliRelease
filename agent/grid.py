from __future__ import absolute_import, division, print_function

import numpy as np

def raster(location, dx):
	return map(int, np.floor(location / dx))

class Line(object):
	def __init__(self, begin, end):
		self.begin = begin
		self.end = end

	def render(self, dx):
		begin = raster(self.begin, dx)[1]
		end = raster(self.end, dx)[1] + 1
		interval = end - begin

		def interpolate(y):
			ratio = float(y - begin) / interval
			return self.begin * (1 - ratio) + self.end * ratio

		return [
			raster(interpolate(y), dx)
			for y in xrange(begin, end + 1)]

class Chain(object):
	def __init__(self, points):
		self.points = points

	def render(self, dx):
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

class Rectangle(object):
	def __init__(self, dimension, location, orientation):
		self.dimension = np.array(dimension)
		self.location = np.array(location)
		self.orientation = orientation

	def corners(self):
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
		perimeter = [left, right]

		indices = np.concatenate([
			map(lambda x: [x, left[1]],
				range(left[0], right[0] + 1))
			for left, right
			in zip(left, right)])

		return set(map(tuple, indices))

class Grid(object):
	def __init__(self, bounds, dx):
		"""
		Create a grid out of the given bounds
		"""

		self.bounds = np.array(bounds)
		self.dx = dx
		self.dimension = self.bounds / dx
		self.grid = np.zeros(map(int, self.dimension))

	def impress(self, indices):
		for index in indices:
			self.grid[index[0]][index[1]] += 1

	def check_in_bounds(self, indices):
		for (x,y) in indices:
			if x >= self.grid.shape[0] or x < 0 or y >= self.grid.shape[0] or y < 0:
				return False
		return True

	def overlap(self):
		total_overlap = np.sum(np.sum(self.grid[self.grid>1]))

		return total_overlap

	def get_forces(self, midpoint, indices):
		location = tuple([int(i / self.dx) for i in midpoint])

		force_vectors = set()
		total_force = [0.0, 0.0]

		for index in indices:
			# find collisions, get force vector for each.
			if self.grid[index] > 1:
				force_vector = tuple([location[0] - index[0], location[1] - index[1]])
				force_vectors.update([force_vector])

		for force in force_vectors:
			total_force[0] += force[0]
			total_force[1] += force[1]

		return total_force