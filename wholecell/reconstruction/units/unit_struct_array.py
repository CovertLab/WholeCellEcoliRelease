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


# TODO: Write test!
class UnitStructArray(object):
	"""UnitStructArray"""

	def __init__(self, struct_array, units):
		self.struct_array = struct_array
		self.units = units

	def field(self, fieldname):
		if self.units[fieldname] == None:
			return self.struct_array[fieldname]
		else:
			return Q_(self.struct_array[fieldname], self.units[fieldname])

	def fieldIs(self, fieldname, new_value, new_units):
		self.structArray[fieldname] = new_value
		self.units[fieldname] = new_units

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
		if type(value) == Q_:
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