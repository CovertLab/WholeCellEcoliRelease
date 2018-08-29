import re
import unum
import types
import numbers
import functools
import collections

import numpy as np
import sympy
from sympy.matrices import dense
import Bio.Seq

import wholecell.utils.unit_struct_array

leaf_types = (
	set,
	str,
	unicode,
	np.int_,
	np.intc,
	np.intp,
	np.int8,
	np.bool_,
	np.dtype,
	np.int16,
	np.int32,
	np.int64,
	np.uint8,
	np.uint16,
	np.uint32,
	np.uint64,
	np.float_,
	np.float16,
	np.float32,
	np.float64,
	np.ndarray,
	np.complex_,
	np.complex64,
	np.complex128,
	unum.Unum,
	Bio.Seq.Seq,
	sympy.Basic,
	numbers.Number,
	functools.partial,
	types.FunctionType,
	dense.MutableDenseMatrix,
	wholecell.utils.unit_struct_array.UnitStructArray)

def is_hidden(attr):
	"""
	Predicate to determine if a given string represents an attribute or function we want to ignore.
	In our case this means that it is flanked by double underscores.
	"""
	return re.search(r'^__.*__$', attr)

def is_leaf(value, leaves=leaf_types):
	"""
	Predicate to determine if we have reached the end of how deep we want to traverse through 
	the object tree. In this case this is true if the object is `callable()` or if it is an 
	instance of one of our leaf types.
	"""
	return callable(value) or isinstance(value, leaves)

def object_tree(obj, path='', debug=None):
	"""
	Given an object, exhaustively traverse down all attributes it contains until leaves are
	reached and convert everything found into a dictionary.

	The resulting dictionary will mirror the structure of the original object, but instead of 
	attributes with values it will be a dictionary where the keys are the attribute names. 
	The type of the dictionarified object will be encoded under the key `!type` which is assumed
	to not be in conflict with any other attributes. The result should aid in serialization and
	deserialization of the object and is intended to be a translation of a pickled object.

	Args:
		obj (object): The object to inspect. 
		path (optional str): The root path of this object tree. This will be built upon for 
			each child of the current object found and reported in a value is provided for `debug`.
		debug (optional str): If provided, prints paths of the attributes encountered. If the
			value is 'ALL', it will print every path. If the value is 'CALLABLE', it will only 
			print methods and functions it finds.
	"""

	if debug == 'ALL':
		print(path)

	if is_leaf(obj):
		if callable(obj) and (debug == 'CALLABLE'):
			print(path)
		return obj
	elif isinstance(obj, collections.Mapping):
		return {key: object_tree(obj[key], "{}['{}']".format(path, key), debug) for key in obj.viewkeys()}
	elif isinstance(obj, collections.Sequence):
		return [object_tree(subobj, "{}[{}]".format(path, index), debug) for index, subobj in enumerate(obj)]
	else:
		attrs = dir(obj)
		tree = {attr: object_tree(getattr(obj, attr), "{}.{}".format(path, attr), debug)
				for attr in attrs
				if not is_hidden(attr)}
		tree['!type'] = type(obj)

		return tree
