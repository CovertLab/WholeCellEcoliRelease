from __future__ import absolute_import, division, print_function

import unum
import types
import numbers
import functools
import collections

import sympy
from sympy.matrices import dense
import Bio.Seq

import wholecell.utils.unit_struct_array


leaf_types = (
	unum.Unum,
	Bio.Seq.Seq,
	sympy.Basic,
	numbers.Number,
	functools.partial,
	types.FunctionType,
	dense.MutableDenseMatrix,
	wholecell.utils.unit_struct_array.UnitStructArray)

def has_python_vars(obj):
	"""
	Returns true if the given object has any Python instance variables, that is
	ordinary fields or compact slots. If not, it's presumably a built-in type
	or extension type implemented entirely in C and Cython.
	"""
	return hasattr(obj, '__dict__') or hasattr(obj, '__slots__')

def all_vars(obj):
	"""
	Returns a dict of all the object's instance variables stored in ordinary
	fields and in compact slots. This expands on the built-in function `vars()`.
	"""
	attrs = getattr(obj, '__dict__', {})
	attrs.update({key: getattr(obj, key) for key in getattr(obj, '__slots__', ())})
	return attrs

def is_leaf(value, leaves=leaf_types):
	"""
	Predicate to determine if we have reached the end of how deep we want to traverse through 
	the object tree.
	"""
	if isinstance(value, (collections.Mapping, collections.Sequence)):
		return isinstance(value, basestring)
	return (callable(value)                 # it's callable
			or isinstance(value, leaves)    # it's an instance of a declared leaf type
			or not has_python_vars(value))  # an object without Python instance variables

def object_tree(obj, path='', debug=None):
	"""
	Diagnostic tool to inspect a complex data structure.

	Given an object, exhaustively traverse down all attributes it contains until leaves are
	reached, and convert everything found into a dictionary or a list.

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
			print('{}: {}'.format(path, obj))
		return obj
	elif isinstance(obj, collections.Mapping):
		return {key: object_tree(value, "{}['{}']".format(path, key), debug)
			for (key, value) in obj.iteritems()}
	elif isinstance(obj, collections.Sequence):
		return [object_tree(subobj, "{}[{}]".format(path, index), debug) for index, subobj in enumerate(obj)]
	else:
		attrs = all_vars(obj)
		tree = {key: object_tree(value, "{}.{}".format(path, key), debug)
				for (key, value) in attrs.iteritems()}
		tree['!type'] = type(obj)

		return tree
