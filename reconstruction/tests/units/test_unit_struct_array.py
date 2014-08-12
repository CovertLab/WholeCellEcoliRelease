'''
test_unit_struct_array.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@self.data: Created 4/18/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.reconstruction.units.unit_struct_array import UnitStructArray
from wholecell.reconstruction.units.unit_registration import Q_


class Test_BulkObjectsContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.data = np.zeros(3, dtype = [('name', 'a10'), ('age', 'int')])
		self.data['name'] = ['nick','derek','john']
		self.data['age'] = [27, 25, 24]
		self.units = {'name' : None, 'age' : 'years'}
		self.usa = UnitStructArray(self.data, self.units)

	def tearDown(self):
		pass

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_init(self):
		with self.assertRaises(Exception) as context:
			UnitStructArray(1., {'hello' : 'goodbye'})
		self.assertEqual(context.exception.message, 'UnitStructArray must be initalized with a numpy array!\n')

		with self.assertRaises(Exception) as context:
			UnitStructArray(self.data, 'foo')
		self.assertEqual(context.exception.message, 'UnitStructArray must be initalized with a dict storing units!\n')

		with self.assertRaises(Exception) as context:
			self.units['hi'] = 'bye'
			UnitStructArray(self.data, self.units)
		self.assertEqual(context.exception.message, 'Struct array fields do not match unit fields!\n')

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_field(self):
		self.assertEqual(
			self.usa._field('name').tolist(),
			['nick','derek','john']
			)

		self.assertTrue(
			type(self.usa._field('age')) == Q_
			)

		self.assertEqual(
			self.usa._field('age').magnitude.tolist(),
			[27, 25, 24]
			)

		self.assertEqual(
			self.usa._field('age').units,
			Q_(0, 'years').units
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_fullArray(self):
		self.assertEqual(
			self.usa.fullArray().tolist(),
			self.data.tolist()
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_fullUnits(self):
		self.assertEqual(
			self.usa.fullUnits(),
			self.units
			)
	
	@noseAttrib.attr('smalltest','unitstructarray')
	def test_getItem_slice(self):
		self.assertEqual(
			self.usa[:2],
			UnitStructArray(self.data[:2], self.units)
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_getItem_indicies(self):
		index = [0,2]

		self.assertEqual(
			self.usa[index],
			UnitStructArray(self.data[index], self.units)
			)

		index = [True, False, True]

		self.assertEqual(
			self.usa[index],
			UnitStructArray(self.data[index], self.units)
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_getItem_field(self):
		self.assertEqual(
			self.usa['name'].tolist(),
			['nick','derek','john']
			)

		self.assertTrue(
			type(self.usa['age']) == Q_
			)

		self.assertEqual(
			self.usa['age'].magnitude.tolist(),
			[27, 25, 24]
			)

		self.assertEqual(
			self.usa['age'].units,
			Q_(0, 'years').units
			)

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_setItem_pintQuantity(self):
		self.usa['age'] = Q_([1,2,3],'year')

		self.assertTrue(
			type(self.usa['age']) == Q_
			)

		self.assertEqual(
			self.usa['age'].magnitude.tolist(),
			[1,2,3]
			)

		self.assertEqual(
			self.usa['age'].units,
			Q_(0, 'years').units
			)

		with self.assertRaises(Exception) as context:
			self.usa['age'] = Q_([1,2,3],'gram')
		self.assertEqual(context.exception.message, 'Units do not match!\n')

	@noseAttrib.attr('smalltest','unitstructarray')
	def test_setItem_structArrayKey(self):
		self.usa['name'] = ['hi','bye','lie']

		self.assertEqual(
			self.usa['name'].tolist(),
			['hi','bye','lie']
			)

		with self.assertRaises(Exception) as context:
			self.usa['age'] = [1,2,3]
		self.assertEqual(context.exception.message, 'Units do not match! Quantity has units your input does not!\n')
		


	# Combo get and set item
	# @noseAttrib.attr('smalltest','unitstructarray')
	# def test_setItem_quantity(self):
	# 	self.usa['age'][:2] = Q_([1,2],'year')
