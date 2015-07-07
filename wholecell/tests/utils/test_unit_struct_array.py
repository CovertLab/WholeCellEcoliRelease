#!/usr/bin/env python

"""
Test unit_struct_array.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 8/14/2014
"""
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.units import g, mol, fg
import numpy as np

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_unit_struct_array(unittest.TestCase):
	

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.struct_array = np.zeros(3, dtype = [('id','a10'),('mass',np.float64)])
		self.units = {'id' : None, 'mass' : g}
		self.us_array = UnitStructArray(self.struct_array, self.units)

	def tearDown(self):
		pass


	@noseAttrib.attr('smalltest','unitstructarray')
	def test_init(self):
		with self.assertRaises(Exception) as context:
			UnitStructArray(1., {'hello' : 'goodbye'})
		self.assertEqual(context.exception.message, 'UnitStructArray must be initialized with a numpy array!\n')

		with self.assertRaises(Exception) as context:
			UnitStructArray(self.struct_array, 'foo')
		self.assertEqual(context.exception.message, 'UnitStructArray must be initialized with a dict storing units!\n')

		with self.assertRaises(Exception) as context:
			self.units['hi'] = 'bye'
			UnitStructArray(self.struct_array, self.units)
		self.assertEqual(context.exception.message, 'Struct array fields do not match unit fields!\n')

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_field(self):
		self.assertTrue(
			self.us_array['id'].tolist(),
			self.struct_array['id'].tolist()
			)

		self.assertTrue(
			(self.us_array['mass'] == g * self.struct_array['mass']).all()
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_fullArray(self):
		self.assertTrue(
			(self.us_array.fullArray() == self.struct_array).all()
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_fullUnits(self):
		self.assertEqual(
			self.us_array.fullUnits(),
			self.units
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_getItem_slice(self):
		self.assertEqual(
			self.us_array[:1],
			UnitStructArray(self.struct_array[:1], self.units)
			)

	@noseAttrib.attr('smalltest','unitstructarray')
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

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_getItem_singleindex(self):
		self.assertEqual(
			self.us_array[0],
			self.struct_array[0]
			)


	@noseAttrib.attr('smalltest','unitstructarray')
	def test_setItem_quantity_with_units(self):
		self.us_array['mass'] = g * np.array([1.,2.,3.])
		self.assertTrue(
			(self.us_array['mass'] == g * np.array([1.,2.,3.])).all()
			)

		with self.assertRaises(Exception) as context:
			self.us_array['mass'] = fg*np.array([1.,2.,3.])
		self.assertEqual(context.exception.message, 'Units do not match!\n')

		with self.assertRaises(Exception) as context:
			self.us_array['mass'] = mol*np.array([1.,2.,3.])
		self.assertEqual(context.exception.message, 'Units do not match!\n')


	@noseAttrib.attr('smalltest','unitstructarray')
	def test_setItem_quantity_no_units(self):
		self.us_array['id'] = ['nick', 'derek', 'john']

		self.assertTrue(
			(self.us_array['id'] == np.array(['nick', 'derek', 'john'])).all()
			)

		with self.assertRaises(Exception) as context:
			self.us_array['mass'] = [1,2,3]
		self.assertEqual(context.exception.message, 'Units do not match! Quantity has units your input does not!\n')
		