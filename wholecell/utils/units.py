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
	if type(array) == Unum:
		units = array.copy()
		units._value = 1
		return units * np.sum(array.asNumber(), axis, dtype, out, keepdims)
	else:
		return np.sum(array, axis, dtype, out, keepdims)


def dot(a, b, out=None):
	if type(a) != Unum:
		a_units = 1
	else:
		a_units = a.copy()
		a_units._value = 1
		a = a.asNumber()

	if type(b) != Unum:
		b_units = 1
	else:
		b_units = b.copy()
		b_units._value = 1
		b  = b.asNumber()
	
	return a_units * b_units * np.dot(a,b,out)

def transpose(array,axis=None):
	if type(a) != Unum or type(b) != Unum:
		raise Exception, 'Only works on Unum!\n'

	units = array.copy()
	units._value = 1

	return units * np.transpose(array.asNumber(), axis)