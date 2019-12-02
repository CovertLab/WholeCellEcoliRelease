"""
Test unit_struct_array.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 8/14/2014
"""

from __future__ import absolute_import, division, print_function

from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.units import g, mol
import numpy as np

import unittest


class Test_unit_struct_array(unittest.TestCase):

	def setUp(self):
		self.struct_array = np.zeros(3, dtype = [('id','a10'),('mass',np.float64)])
		self.units = {'id' : None, 'mass' : g}
		self.us_array = UnitStructArray(self.struct_array, self.units)

	def test_init(self):
		with self.assertRaisesRegexp(
			Exception,
			'^UnitStructArray must be initialized with a numpy array!\n$',
		):
			UnitStructArray(1., {'hello': 'goodbye'})

		with self.assertRaisesRegexp(
			Exception,
			'^UnitStructArray must be initialized with a dict storing '
			'units!\n$',
		):
			UnitStructArray(self.struct_array, 'foo')

		with self.assertRaisesRegexp(
			Exception,
			'Struct array fields do not match unit fields!\n',
		):
			self.units['hi'] = 'bye'
			UnitStructArray(self.struct_array, self.units)

	def test_field(self):
		self.assertEqual(
			self.us_array['id'].tolist(),
			self.struct_array['id'].tolist()
			)

		self.assertTrue(
			(self.us_array['mass'] == g * self.struct_array['mass']).all()
			)

	def test_fullArray(self):
		self.assertTrue(
			(self.us_array.fullArray() == self.struct_array).all()
			)

	def test_fullUnits(self):
		self.assertEqual(
			self.us_array.fullUnits(),
			self.units
			)

	def test_getItem_slice(self):
		self.assertEqual(
			self.us_array[:1],
			UnitStructArray(self.struct_array[:1], self.units)
			)

	def test_getItem_indicies(self):
		index = [0,2]

		self.assertEqual(
			self.us_array[index],
			UnitStructArray(self.struct_array[index], self.units)
			)

		index = [True, False, True]

		self.assertEqual(
			self.us_array[index],
			UnitStructArray(self.struct_array[index], self.units)
			)

	def test_getItem_singleindex(self):
		self.assertEqual(
			self.us_array[0],
			self.struct_array[0]
			)


	def test_setItem_quantity_with_units(self):
		self.us_array['mass'] = g * np.array([1., 2., 3.])
		self.assertTrue(
			(self.us_array['mass'] == g * np.array([1., 2., 3.])).all()
		)

		with self.assertRaisesRegexp(
			Exception,
			'Units do not match!\n',
		):
			self.us_array['mass'] = mol * np.array([1., 2., 3.])


	def test_setItem_quantity_no_units(self):
		self.us_array['id'] = ['nick', 'derek', 'john']

		self.assertTrue(
			(self.us_array['id'] == np.array(['nick', 'derek', 'john'])).all()
		)

		with self.assertRaisesRegexp(
			Exception,
			'Units do not match! Quantity has units your input does not!\n',
		):
			self.us_array['mass'] = [1, 2, 3]
