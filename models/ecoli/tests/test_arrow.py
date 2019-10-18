from __future__ import absolute_import, division, print_function

import numpy as np
from pprint import pprint
import unittest

from arrow import StochasticSystem


class TestArrow(unittest.TestCase):
	"""Test linkage to the stochastic-arrow library."""

	def test_arrow(self):
		stoichiometric_matrix = np.array([
			[1, 1, -1, 0],
			[-2, 0, 0, 1],
			[-1, -1, 1, 0]], np.int64)

		rates = np.array([3, 1, 1]) * 0.01

		arrow = StochasticSystem(stoichiometric_matrix, rates)
		result = arrow.evolve(1.0, np.array([50, 20, 30, 40], np.int64))
		pprint(result)

		self.assertEqual(stoichiometric_matrix.shape[0], arrow.obsidian.reactions_count())
		self.assertEqual(stoichiometric_matrix.shape[1], arrow.obsidian.substrates_count())

		self.assertEqual(13, result['steps'])
		self.assertEqual(13, len(result['events']))
		self.assertEqual(2, result['events'][0])


if __name__ == '__main__':
	unittest.main()
