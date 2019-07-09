from __future__ import absolute_import, division, print_function

import unum
import numpy as np
import re
import types
import numbers
import functools
import collections

import sympy
from sympy.matrices import dense
import Bio.Seq

import wholecell.utils.unit_struct_array

NULP = 0  # float comparison tolerance, in Number of Units in the Last Place

leaf_types = (
	unum.Unum,
	Bio.Seq.Seq,
	sympy.Basic,
	numbers.Number,
	functools.partial,
	types.FunctionType,
	dense.MutableDenseMatrix,
	wholecell.utils.unit_struct_array.UnitStructArray)


WHITESPACE = re.compile(r'\s+')

class Repr(object):
	'''A Repr has the given repr() string without quotes and != any other value.'''
	def __init__(self, repr_):
		self.repr_ = repr_

	def __repr__(self):
		return self.repr_


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
	If the object implements the pickling method `__getstate__`, call that
	instead to get its defining state.
	"""
	if hasattr(obj, '__getstate__'):
		return obj.__getstate__()

	attrs = getattr(obj, '__dict__', {})
	attrs.update({key: getattr(obj, key) for key in getattr(obj, '__slots__', ())})
	return attrs

def is_leaf(value, leaves=leaf_types):
	"""
	Predicate to determine if we have reached the end of how deep we want to traverse
	through the object tree.
	"""
	if isinstance(value, (collections.Mapping, collections.Sequence)):
		return isinstance(value, basestring)
	return (callable(value)                 # it's callable
			or isinstance(value, leaves)    # it's an instance of a declared leaf type
			or not has_python_vars(value))  # an object without Python instance variables

def object_tree(obj, path='', debug=None):
	"""
	Diagnostic tool to inspect a complex data structure.

	Given an object, exhaustively traverse down all attributes it contains until leaves
	are	reached, and convert everything found into a dictionary or a list. The resulting
	dictionary will mirror the structure of the original object, but instead of
	attributes with values it will be a dictionary where the keys are the attribute
	names. The type of the dictionarified object will be encoded under the key `!type`
	which is assumed to not be in conflict with any other attributes. The result should
	aid in serialization and deserialization of the object and is intended to be a
	translation of a pickled object.

	Args:
		obj (object): The object to inspect. 
		path (optional str): The root path of this object tree. This will be built upon
	        for each child of the current object found and reported in a value is
	        provided for `debug`.
		debug (optional str): If provided, prints paths of the attributes encountered.
	        If the value is 'ALL', it will print every path. If the value is 'CALLABLE',
	        it will only print methods and functions it finds.
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

def diff_trees(a, b):
	"""
	Find the differences between two trees or leaf nodes a and b. Return a
	falsely value if the inputs match OR a truthy value that explains or
	summarizes their differences, where each point in the tree where the inputs
	differ will be a tuple (a's value, b's value, optional description).

	Floating point numbers are compared with the tolerance set by the constant
	NULP (Number of Units in the Last Place), allowing for NaN and infinite
	values. (Adjust the tolerance level NULP if needed.)

	This operation is symmetrical.
	"""

	# if they aren't they same type, they are clearly different. Also this lets us
	# safely assume throughout the rest of the function that a and b are the same type
	if type(a) != type(b):
		return elide(a, max_len=400), elide(b, max_len=400)

	# if they are floats, handle various kinds of values
	elif isinstance(a, float):
		return compare_floats(a, b)

	# if they are numpy arrays, compare them using a numpy testing function
	elif isinstance(a, np.ndarray):
		return compare_ndarrays(a, b)

	# if they are Unums compare their contents with matching units
	elif isinstance(a, unum.Unum):
		a0, b0 = a.matchUnits(b)
		return diff_trees(a0.asNumber(), b0.asNumber())

	# if they are leafs (including strings) use python equality comparison
	elif is_leaf(a):
		if a != b:
			return elide(a), elide(b)

	# if they are dictionaries then diff the value under each key
	elif isinstance(a, collections.Mapping):
		diff = {}
		na = Repr('--')
		nb = Repr('--')
		for key in set(a.keys()) | set(b.keys()):
			subdiff = diff_trees(a.get(key, na), b.get(key, nb))
			if subdiff:
				diff[key] = subdiff
		return diff

	# if they are sequences then compare each element at each index
	elif isinstance(a, collections.Sequence):
		if len(a) > len(b):
			b = list(b) + (len(a) - len(b)) * [Repr('--')]
		elif len(b) > len(a):
			a = list(a) + (len(b) - len(a)) * [Repr('--')]

		diff = []
		for index in xrange(len(a)):
			subdiff = diff_trees(a[index], b[index])
			if subdiff:
				diff.append(subdiff)
		return diff

	# this should never happen
	else:
		print('value not considered by `diff_trees`: {} {}'.format(a, b))

def elide(value, max_len=200):
	'''Return a value with the same repr but elided if it'd be longer than max.'''
	repr_ = repr(value)
	if len(repr_) > max_len:
		return Repr(repr_[:max_len] + '...')
	return value

def simplify_error_message(message):
	return elide(Repr(WHITESPACE.sub(' ', message).strip()))

def compare_floats(f1, f2):
	'''Compare two floats, allowing some tolerance, NaN, and Inf values. This
	considers all types of NaN to match.
	Return 0.0 (which is falsey) if they match, else (f1, f2).
	'''
	if f1 == f2 or np.isnan(f1) and np.isnan(f2):
		return 0.0
	try:
		np.testing.assert_array_almost_equal_nulp(f1, f2, nulp=NULP)
		return 0.0
	except AssertionError:
		# FWIW, the error.message tells the NULP difference.
		return f1, f2

def compare_ndarrays(array1, array2):
	'''Compare two ndarrays, checking the shape and all elements, allowing for
	NaN values and non-numeric values. Return () if they match, else a tuple of
	diff info or just a diff description.

	TODO(jerry): Allow tolerance for float elements of structured arrays and
	  handle NaN and Inf values.
	'''
	if issubclass(array1.dtype.type, np.floating):
		try:
			# This handles float tolerance but not NaN and Inf.
			np.testing.assert_array_almost_equal_nulp(array1, array2, nulp=NULP)
			return ()
		except AssertionError as e:
			# return elide(array1), elide(array2), simplify_error_message(e.message)
			pass  # try again, below

	try:
		# This handles non-float dtypes, also NaN and Inf, but no tolerance.
		np.testing.assert_array_equal(array1, array2)
		return ()
	except AssertionError as e:
		return simplify_error_message(e.message)
