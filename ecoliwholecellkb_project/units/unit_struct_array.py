#!/usr/bin/env python

"""
UnitStructArray

Wraps Numpy struct arrays using Pint units. Will assure that correct units are
being used while loading constants.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/31/2014
"""

from pint import UnitRegistry
UREG = UnitRegistry()
Q_ = UREG.Quantity

class UnitStructArray(object):
	"""UnitStructArray"""

	def __init__(self, struct_array, units):
		self.struct_array = struct_array
		self.units = units

	def field(self, fieldname):
		return Q_(self.struct_array[fieldname], self.units[fieldname])

	def fieldIs(self, fieldname, new_value, new_units):
		self.structArray[fieldname] = new_value
		self.units[fieldname] = new_units

	def fullArray(self):
		return self.struct_array

	def fullUnits(self):
		return self.units