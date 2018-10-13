from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
import numpy.testing as npt
import nose.plugins.attrib as noseAttrib

from agent.grid import Line, Chain, Rectangle, Grid

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


	@noseAttrib.attr('grid')
	def test_line(self):
		print('line')
		print(self.line.render(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid')
	def test_chain(self):
		print('chain')
		print(self.chain.render(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid')
	def test_rectangle(self):
		print('rectangle')
		print(self.rectangle.render(0.1))

		npt.assert_equal(1, 1)

	@noseAttrib.attr('grid')
	def test_grid(self):
		self.grid.impress(self.rectangle)

		print('grid')
		unique, counts = np.unique(np.concatenate(self.grid.grid), return_counts=True)
		print(dict(zip(unique, counts)))

		npt.assert_equal(1, 1)
