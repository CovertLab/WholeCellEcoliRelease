"""
Miscellaneous data structure utilities.
"""

from __future__ import absolute_import, division, print_function

import os
from typing import Any, Dict, Iterable, Mapping

from six import string_types


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

def expand_keyed_env_vars(data):
	# type: (Any) -> Any
	"""If data is a dict, return a new dict where each _keyed environment
	'$VARIABLE' value is expanded into a non-underscore entry, e.g.:

	{'_foo': '$HOME'} --> {'foo': '/Users/franklin', '_foo': '$HOME'}

	To trigger expansion, an entry must have a key string starting with '_' and
	a value string starting with '$' that names an environment variable.
	"""
	if not isinstance(data, dict):
		return data

	copy = data.copy()
	for key, val in data.items():
		if startswith(key, '_') and startswith(val, '$'):
			expanded = os.environ.get(val[1:])
			if expanded:
				copy[key[1:]] = expanded

	return copy

def select_keys(mapping, keys, **kwargs):
	# type: (Mapping[str, Any], Iterable[str], **Any) -> Dict[str, Any]
	"""Return a dict of the selected keys from mapping plus the kwargs."""
	result = {key: mapping[key] for key in keys}
	result.update(**kwargs)
	return result

def startswith(obj, prefix):
	# type: (Any, str) -> bool
	"""Return True if obj is a string that starts with the given prefix."""
	return isinstance(obj, string_types) and obj.startswith(prefix)
