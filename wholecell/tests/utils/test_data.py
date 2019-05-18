"""Test the data utility."""

from __future__ import absolute_import, division, print_function

import nose.plugins.attrib as noseAttrib
import nose.tools
import unittest

from wholecell.utils import data


@noseAttrib.attr('smalltest', 'data')
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

		with nose.tools.assert_raises(KeyError):
			data.dissoc_strict(source, (100, 200))

	def test_select_keys(self):
		source = {i: i * 10 for i in range(5)}
		d1 = data.select_keys(source, [1, 4])
		self.assertEqual(d1, {1: 10, 4: 40})

		with nose.tools.assert_raises(KeyError):
			data.select_keys(source, (2, 100))


if __name__ == '__main__':
	unittest.main()
