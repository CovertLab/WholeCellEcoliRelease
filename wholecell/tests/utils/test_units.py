from __future__ import absolute_import, division, print_function

import unittest
import nose.plugins.attrib as noseAttrib

from wholecell.utils import units


@noseAttrib.attr('smalltest', 'units')
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

	def test_slash_rdiv(self):
		"""Test __rtruediv__ by virtue of __future__.division ."""
		x = 1
		y = 100 * units.s
		quotient = (x / y).asNumber()
		self.assertEqual(0.01, quotient)

	# TODO(jerry): Test the array functions.
