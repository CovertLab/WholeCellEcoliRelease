from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
import numpy.testing as npt
import nose.plugins.attrib as noseAttrib

from environment.collision.grid import Line, Chain, Rectangle, Grid

class TestGrid(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		a = np.array([1.0, 2.0])
		b = np.array([4.0, 4.0])
		c = np.array([5.0, 7.0])

		self.line = Line(a,b)
		self.chain = Chain([a, b, c])
		self.rectangle = Rectangle(a, b, np.pi / 5)
		self.grid = Grid(np.array([10.0, 10.0]), 0.1)


	def tearDown(self):
		pass


	@noseAttrib.attr('grid', 'shape')
	def test_line(self):
		indexes = self.line.indexes(0.1)

		print('line')
		print(indexes)

		npt.assert_equal(len(indexes), 22)

	@noseAttrib.attr('grid', 'shape')
	def test_chain(self):
		indexes = self.chain.indexes(0.1)

		print('chain')
		print(indexes)

		npt.assert_equal(len(indexes), 53)

	@noseAttrib.attr('grid', 'shape')
	def test_rectangle(self):
		indexes = self.rectangle.indexes(0.1)

		print('rectangle')
		print(indexes)

		npt.assert_equal(len(indexes), 229)

	@noseAttrib.attr('grid', 'shape')
	def test_grid(self):
		self.grid.impress(self.rectangle)

		unique, counts = np.unique(np.concatenate(self.grid.grid), return_counts=True)
		outcome = dict(zip(unique, counts))

		print('grid')
		print(outcome)

		npt.assert_equal(outcome[0], 229)
