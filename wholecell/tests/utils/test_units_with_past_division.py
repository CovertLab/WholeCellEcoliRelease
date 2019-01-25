from __future__ import absolute_import, print_function  # <=== NOT division!

import sys
import unittest
import nose.plugins.attrib as noseAttrib

from wholecell.utils import units

PY2 = sys.version_info.major < 3


@noseAttrib.attr('smalltest', 'units')
class Test_units_with_past_division(unittest.TestCase):

	def test_slash_div(self):
		"""Test Python 2 `/` *WITHOUT* __future__.division."""
		x = 1 * units.s
		y = 100 * units.s
		quotient = (x / y).asNumber()
		expected = 0 if PY2 else 0.01
		self.assertEqual(expected, quotient)

	def test_slash_rdiv(self):
		"""Test Python 2 __rdiv__ `/` *WITHOUT* __future__.division."""
		x = 1
		y = 100 * units.s
		quotient = (x / y).asNumber()
		expected = 0 if PY2 else 0.01
		self.assertEqual(expected, quotient)
