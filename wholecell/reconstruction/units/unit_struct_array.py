#!/usr/bin/env python

"""
UnitStructArray

Wraps Numpy struct arrays using Pint units. Will assure that correct units are
being used while loading constants.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/31/2014
"""

from unit_registration import Q_
import numpy as np

class UnitStructArray(object):
	"""UnitStructArray"""

	def __init__(self, struct_array, units):
		self._validate(struct_array, units)

		self.struct_array = struct_array
		self.units = units

	def _validate(self, struct_array, units):
		s = ''
		if type(struct_array) != np.ndarray:
			s += 'UnitStructArray must be initalized with a numpy array!\n'
		elif type(units) != dict:
			s += 'UnitStructArray must be initalized with a dict storing units!\n'
		elif set([x[0] for x in struct_array.dtype.descr]) != set(units.keys()):
			s += 'Struct array fields do not match unit fields!\n'
		if len(s):
			raise Exception, s

	def _field(self, fieldname):
		if self.units[fieldname] == None:
			return self.struct_array[fieldname]
		else:
			return Q_(self.struct_array[fieldname], self.units[fieldname])

	def fullArray(self):
		return self.struct_array

	def fullUnits(self):
		return self.units

	def __getitem__(self, key):
		if type(key) == slice:
			return UnitStructArray(self.struct_array[key], self.units)
		elif type(key) == np.ndarray or type(key) == list:
			return UnitStructArray(self.struct_array[key], self.units)
		else:
			return self.field(key)

	def __setitem__(self, key, value):
		if type(value) == pint.unit.Quantity:
			if self.units[key] != value.units:
				raise Exception, 'Units do not match!\n'
			self.structArray[key] = value.magnitude
			self.units[key] = value.units
		else:
			if self.units[key] != None:
				raise Exception, 'Units do not match! Quantity has units your input does not!\n'
			self.structArray[key] = value
			self.units[key] = None

	def __len__(self):
		return len(self.struct_array)

	def __repr__(self):
		return 'STRUCTURED ARRAY:\n{}\nUNITS:\n{}'.format(self.struct_array.__repr__(), self.units)