#!/usr/bin/env python

"""
Units

Defines/registers custom units for Pint

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/14/2014
"""

import scipy.constants
import numpy as np
from unum.units import *
from unum import Unum

count = Unum.unit('count',mol/(scipy.constants.Avogadro))
nt = Unum.unit('nucleotide', count)
aa = Unum.unit('amino_acid', count)


def sum(array, axis = None, dtype=None, out=None, keepdims=False):
	if not isinstance(array,Unum):
		raise Exception("Only works on Unum!")

	units = getUnit(array)
	return units * np.sum(array.asNumber(), axis, dtype, out, keepdims)		

def dot(a, b, out=None):
	if not isinstance(a, Unum):
		a_units = 1
	else:
		a_units = getUnit(a)
		a = a.asNumber()

	if not isinstance(b,Unum):
		b_units = 1
	else:
		b_units = getUnit(b)
		b  = b.asNumber()
	
	return a_units * b_units * np.dot(a,b,out)

def floor(x):
	if not hasUnit(x):
		raise Exception('Only works on Unum!')
	x_unit = getUnit(x)
	x = x.asNumber()
	return x_unit * np.floor(x)

def transpose(array,axis=None):
	if not isinstance(a,Unum):
		raise Exception('Only works on Unum!')
	if not isinstance(b,Unum):
		raise Exception('Only works on Unum!')

	units = getUnit(array)

	return units * np.transpose(array.asNumber(), axis)

def hstack(tup):
	unit = getUnit(tup[0])
	value = []
	for array in tup:
		if not isinstance(array,Unum):
			raise Exception('Only works on Unum!')
		else:
			array.normalize()
			value.append(array.matchUnits(unit)[0].asNumber())
	value = tuple(value)
	return unit * np.hstack(value)

def getUnit(value):
	if not hasUnit(value):
		raise Exception("Only works on Unum!")

	value.normalize()
	value_units = value.copy()
	value_units._value = 1
	return value_units

def hasUnit(value):
	if isinstance(value, Unum):
		return True
	else:
		return False

def convertNoUnitToNumber(value):
	if not hasUnit(value):
		raise Exception("Only works on Unum!")

	value.normalize()
	value.checkNoUnit()
	return value.asNumber()