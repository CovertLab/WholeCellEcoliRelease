"""
Miscellaneous data structure utilities.
"""

from __future__ import absolute_import, division, print_function

from typing import Iterable, Mapping


def dissoc(mapping, keys):
	# type: (Mapping, Iterable) -> dict
	"""Dissociate: Return a new dict like `mapping` without the given `keys`.
	See also `dissoc_strict()`.
	"""
	result = dict(mapping)

	for key in keys:
		result.pop(key, None)
	return result

def dissoc_strict(mapping, keys):
	# type: (Mapping, Iterable) -> dict
	"""Dissociate: Return a new dict like `mapping` without the given `keys`.
	Raises a KeyError if any of these keys are not in the given mapping.
	"""
	result = dict(mapping)

	for key in keys:
		del result[key]
	return result
