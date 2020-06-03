"""
UnitStructArray

Wraps Numpy struct arrays using Pint units. Will assure that correct units are
being used while loading constants.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/31/2014
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from typing import Dict, Optional
import unum

from wholecell.utils import units as units_pkg

class UnitStructArray(object):
	"""Wraps Numpy structured arrays using Pint units. Will assure that correct
	units are being used while loading constants.
	"""

	def __init__(self, struct_array, units):
		# type: (np.ndarray, Dict[str, Optional[str]]) -> None
		self._validate(struct_array, units)

		self.struct_array = struct_array
		self.units = units

	def _validate(self, struct_array, units):
		s = ''
		if not isinstance(struct_array, np.ndarray):
			s += 'UnitStructArray must be initialized with a numpy array!\n'
		elif not isinstance(units, dict):
			s += 'UnitStructArray must be initialized with a dict storing units!\n'
		elif set([x[0] for x in struct_array.dtype.descr]) != set(units.keys()):
			s += 'Struct array fields do not match unit fields!\n'
		if len(s):
			raise Exception(s)

	def _field(self, fieldname):
		if not units_pkg.hasUnit(self.units[fieldname]):
			if self.units[fieldname] is None:
				return self.struct_array[fieldname]
			else:
				raise Exception('Field has incorrect units or unitless designation!\n')
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
		if units_pkg.hasUnit(value):
			try:
				self.units[key].matchUnits(value)
			except unum.IncompatibleUnitsError:
				raise Exception('Units do not match!\n')

			self.struct_array[key] = value.asNumber()
			self.units[key] = units_pkg.getUnit(value)

		elif type(value) == list or type(value) == np.ndarray:
			if units_pkg.hasUnit(self.units[key]):
				raise Exception('Units do not match! Quantity has units your input does not!\n')
			self.struct_array[key] = value
			self.units[key] = None

		else:
			raise Exception("Can't assign data-type other than unum datatype or list/numpy array!\n")

	def __len__(self):
		return len(self.struct_array)

	def __repr__(self):
		return 'STRUCTURED ARRAY:\n{}\nUNITS:\n{}'.format(self.struct_array.__repr__(), self.units)

	def __eq__(self, other):
		if type(other) != type(self):
			return False
		elif self.struct_array.dtype != other.struct_array.dtype:
			return False
		elif not all(self.struct_array == other.struct_array):
			return False
		elif self.units != other.units:
			return False
		else:
			return True

	def __ne__(self, other):
		return not self.__eq__(other)
