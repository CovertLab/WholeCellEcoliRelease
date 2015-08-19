#!/usr/bin/env python

"""
UnitStructArray

Wraps Numpy struct arrays using Pint units. Will assure that correct units are
being used while loading constants.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/31/2014
"""

import numpy as np
import unum


# TODO: Write test!
class UnitStructArray(object):
	"""UnitStructArray"""

	def __init__(self, struct_array, units):
		self._validate(struct_array, units)

		self.struct_array = struct_array
		self.units = units

	def _validate(self, struct_array, units):
		s = ''
		if type(struct_array) != np.ndarray:
			s += 'UnitStructArray must be initialized with a numpy array!\n'
		elif type(units) != dict:
			s += 'UnitStructArray must be initialized with a dict storing units!\n'
		elif set([x[0] for x in struct_array.dtype.descr]) != set(units.keys()):
			s += 'Struct array fields do not match unit fields!\n'
		if len(s):
			raise Exception, s

	def _field(self, fieldname):
		if type(self.units[fieldname]) != unum.Unum:
			if self.units[fieldname] == None:
				return self.struct_array[fieldname]
			else:
				raise Exception, 'Field has incorrect units or unitless designation!\n'
		else:
			return self.units[fieldname] * self.struct_array[fieldname]

	def fullArray(self):
		return self.struct_array

	def fullUnits(self):
		return self.units

	def __getitem__(self, key):
		if type(key) == slice:
			return UnitStructArray(self.struct_array[key], self.units)
		elif type(key) == np.ndarray or type(key) == list:
			return UnitStructArray(self.struct_array[key], self.units)
		elif type(key) == int:
			return self.struct_array[key]
		else:
			return self._field(key)

	def __setitem__(self, key, value):
		if type(value) == unum.Unum:
			if self.units[key].strUnit() != value.strUnit():
				raise Exception, 'Units do not match!\n'
			self.struct_array[key] = value.asNumber()
			# This is a bit of a hack but I couldn't figure out a
			# method to get just the units object out of a Unum
			value_units = value.copy()
			value_units._value = 1
			self.units[key] = value_units
		elif type(value) == list or type(value) == np.ndarray:
			if type(self.units[key]) == unum.Unum:
				raise Exception, 'Units do not match! Quantity has units your input does not!\n'
			self.struct_array[key] = value
			self.units[key] = None
		else:
			raise Exception, 'Cant assign data-type other than unum datatype or list/numpy array!\n'

	def __len__(self):
		return len(self.struct_array)

	def __repr__(self):
		return 'STRUCTURED ARRAY:\n{}\nUNITS:\n{}'.format(self.struct_array.__repr__(), self.units)

	def __eq__(self, other):
		if type(other) != type(self):
			return False
		elif not all(self.struct_array == other.struct_array):
			return False
		elif self.units != other.units:
			return False
		else:
			return True

	def __ne__(self, other):
		return not self.__eq__(other)
