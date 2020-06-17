from __future__ import absolute_import, division, print_function

import unittest
import numpy as np

from wholecell.utils import units


class Test_units(unittest.TestCase):

	def test_basics(self):
		"""Test abs(), floor(), hasUnit(), getUnit()."""
		m3 = 3 * units.m
		m314 = 3.14 * units.m

		self.assertEqual(m3, units.floor(m314))

		self.assertEqual(m314, units.abs(m314))
		self.assertEqual(m314, units.abs(-m314))

		self.assertFalse(units.hasUnit(3))
		self.assertTrue(units.hasUnit(m3))
		self.assertEqual(units.m, units.getUnit(m3))

		with self.assertRaises(Exception):
			units.getUnit(3)

	def test_slash_div(self):
		"""Test that `/` with __future__.division does __truediv__"""
		x = 1 * units.s
		y = 100 * units.s
		quotient = (x / y).asNumber()

		self.assertEqual(0.01, quotient)

	def test_truediv(self):
		"""Test __truediv__."""
		x = 1 * units.s
		y = 100 * units.s
		quotient = x.__truediv__(y).asNumber()

		self.assertEqual(0.01, quotient)

	def test_div_int(self):
		"""Test __div__ with ints."""
		x = 1 * units.s
		y = 100 * units.s
		quotient = x.__div__(y).asNumber()
		self.assertEqual(0, quotient)

	def test_div_float(self):
		"""Test __div__ with float and int."""
		x = 1.0 * units.s
		y = 100 * units.s
		quotient = x.__div__(y).asNumber()
		self.assertEqual(0.01, quotient)

		quotient = (x / y).asNumber()
		self.assertEqual(0.01, quotient)

		y = 100
		quotient = x.__div__(y).asNumber()
		self.assertEqual(0.01, quotient)

		quotient = (x / y).asNumber()
		self.assertEqual(0.01, quotient)

	def test_floordiv(self):
		"""Test __floordiv__ via `//`."""
		x = 150 * units.s
		y = 100
		quotient = (x // y).asNumber()
		self.assertEqual(1, quotient)

		d = np.array([1, 2, 3])
		d2 = 2 // d
		e1 = units.s * (2 // d)
		e2 = units.s * 2 // d
		np.testing.assert_array_equal(d2, e1.asNumber())
		np.testing.assert_array_equal(e1, e2)

	def test_rfloordiv(self):
		"""Test __rfloordiv__ via `//`."""
		x = 150
		y = 100 * units.s
		quotient = (x // y).asNumber()  # __rfloordiv__
		self.assertEqual(1, quotient)

		d = np.array([1, 2, 3])
		e1 = 2 // d
		e2 = units.s * 2 // d           # __floordiv__
		a1 = 2 // (1 // units.s * d)    # __rfloordiv__
		np.testing.assert_array_equal(e1, e2.asNumber())
		np.testing.assert_array_equal(e2, a1)

	def test_slash_rdiv(self):
		"""Test __rtruediv__ by via `/` and __future__.division ."""
		x = 1
		y = 100 * units.s
		quotient = (x / y).asNumber()   # __rtruediv__

		self.assertEqual(0.01, quotient)

		d = np.array([1, 2, 3])
		d2 = 2 / d
		e1 = units.s * (2 / d)          # __truediv__
		e2 = units.s * 2 / d            # __truediv__
		a1 = 2 / (1 / units.s * d)      # __rtruediv__
		a2 = 2 / (units.s * d)          # __rtruediv__

		np.testing.assert_array_equal(e1, e2)
		np.testing.assert_array_equal(e1, a1)
		np.testing.assert_array_equal(d2, a2.asNumber())

	# TODO(jerry): Test the array functions.
