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
		print('line')
		print(self.line.indexes(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid', 'shape')
	def test_chain(self):
		print('chain')
		print(self.chain.indexes(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid', 'shape')
	def test_rectangle(self):
		print('rectangle')
		print(self.rectangle.indexes(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid', 'shape')
	def test_grid(self):
		self.grid.impress(self.rectangle)

		unique, counts = np.unique(np.concatenate(self.grid.grid), return_counts=True)
		outcome = dict(zip(unique, counts))

		print('grid')
		print(outcome)

		npt.assert_equal(outcome[0], 229)
