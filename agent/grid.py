from __future__ import absolute_import, division, print_function

import numpy as np

def raster(location, dx):
	return map(int, np.floor(location / dx))

def within(bounds, point):
	return point[0] >= 0 and point[1] >= 0 and point[0] < bounds[0] and point[1] < bounds[1]

class Shape(object):
	def __init__(self):
		self.renders = {}

	def indices(self, dx):
		if dx not in self.renders:
			self.renders[dx] = self.render(dx)
		return self.renders[dx]

class Line(Shape):
	def __init__(self, begin, end):
		super(Line, self).__init__()

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

class Chain(Shape):
	def __init__(self, points):
		super(Chain, self).__init__()

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

class Rectangle(Shape):
	def __init__(self, dimension, location, orientation):
		super(Rectangle, self).__init__()

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
		self.grid = np.full(map(int, self.dimension), -1)

	def reset(self):
		self.grid.fill(-1)

	def impress(self, shape):
		indices = shape.indices(self.dx)
		for index in indices:
			if within(self.dimension, index):
				self.grid[index[0]][index[1]] += 1

	def check_in_bounds(self, shape):
		for index in shape.indices(self.dx):
			if within(self.grid.shape, index):
				return False
		return True

	def overlap(self):
		total_overlap = np.sum(np.sum(self.grid[self.grid > 0]))
		return total_overlap

	def forces(self, midpoint, shape):
		location = np.array([int(i / self.dx) for i in midpoint])
		total_force = np.array([0.0, 0.0])

		pixels = 0.
		for index in shape.indices(self.dx):
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

		# magnitude = np.linalg.norm(total_force)
		# if magnitude == 0:
		# 	total_force *= 0
		# else:
		# 	total_force /= magnitude
		# total_force *= 15 * (self.dx ** 2)
			
		scaling = 1. / np.sqrt(pixels) if pixels > 0 else 1
		total_force *= scaling * 1 * (self.dx ** 2)

		return total_force
