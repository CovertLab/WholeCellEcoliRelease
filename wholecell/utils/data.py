"""
Miscellaneous data structure utilities.
"""

from __future__ import absolute_import, division, print_function

from typing import Any, Iterable, Mapping


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

def select_keys(mapping, keys, **kwargs):
	# type: (Mapping[str, Any], Iterable[str], **Any) -> Dict[str, Any]
	"""Return a dict of the selected keys from mapping plus the kwargs."""
	result = {key: mapping[key] for key in keys}
	result.update(**kwargs)
	return result
