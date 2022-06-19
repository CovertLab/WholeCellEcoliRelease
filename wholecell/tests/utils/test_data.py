"""Test the data utility."""

from __future__ import absolute_import, division, print_function

import unittest

import pytest

from wholecell.utils import data
from six.moves import range


class Test_data(unittest.TestCase):

	def test_dissoc_and_dissoc_strict(self):
		source = {i: i * 10 for i in range(5)}
		d1 = data.dissoc(source, [1, 3])
		self.assertEqual(d1, {i: i * 10 for i in [0, 2, 4]})

		d2 = data.dissoc_strict(source, [1, 3])
		self.assertEqual(d2, {i: i * 10 for i in [0, 2, 4]})

		d1 = data.dissoc(source, (100, 200))
		self.assertEqual(d1, source)
		self.assertIsNot(d1, source)

		with pytest.raises(KeyError):
			data.dissoc_strict(source, (100, 200))

	def test_select_keys(self):
		source = {'x{}'.format(i): i * 10 for i in range(5)}
		d1 = data.select_keys(source, ['x1', 'x4'])
		self.assertEqual({'x1': 10, 'x4': 40}, d1)

		# Expect absent keys to raise KeyError.
		with pytest.raises(KeyError):
			data.select_keys(source, ('x2', '100'))

		# Test added keys.
		d2 = data.select_keys(source, ('x1', 'x3', 'x4'), x=-1, y=-2, x3=0)
		self.assertEqual({'x1': 10, 'x3': 0, 'x4': 40, 'x': -1, 'y': -2}, d2)


if __name__ == '__main__':
	unittest.main()
