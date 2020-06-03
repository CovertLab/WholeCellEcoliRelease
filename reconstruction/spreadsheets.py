"""
Subclasses of DictWriter and DictReader that parse plaintext as JSON strings,
allowing for basic type parsing and fields that are dictionaries or lists, with
optional units.
"""

from __future__ import absolute_import, division, print_function

import csv
import io
from itertools import ifilterfalse
import json
import re
import numpy as np
from typing import Any, cast, Dict, Sequence, Text

import six

from wholecell.utils import units


CSV_DIALECT = csv.excel_tab


def comment_line(line):
	# type: (str) -> bool
	return line.lstrip().startswith('#')

def read_tsv(filename):
	# type: (str) -> Sequence[Dict[str, Any]]
	"""Read an entire .tsv file using JsonReader/DictReader."""
	# ########################################################################
	# NOTE: Python 3 csv requires opening the file as text 'r' with the right
	# character encoding while Python 2 csv requires opening it as bytes 'rb'
	# then decoding from UTF-8 after csv reads it but before DictReader or at
	# least before json.loads().
	#
	# Several of the .tsv files are in UTF-8. The rest are in the ASCII subset.
	# ########################################################################

	mode = 'rb' if six.PY2 else 'r'
	encoding = None if six.PY2 else 'utf-8'
	with io.open(filename, mode=mode, encoding=encoding) as fh:
		reader = JsonReader(ifilterfalse(comment_line, fh), dialect=CSV_DIALECT)
		return list(reader)


def array_to_list(value):
	if isinstance(value, np.ndarray):
		value = value.tolist()

	return value


class JsonWriter(csv.DictWriter, object):
	# [Python 2 DictWriter is an old-style class so mix in `object` to get a
	# new-style class that supports `super()`.]
	def __init__(self, *args, **kwargs):
		"""Defaults dialect=CSV_DIALECT."""
		kwargs.setdefault('dialect', CSV_DIALECT)
		super(JsonWriter, self).__init__(
			quotechar = "'", quoting = csv.QUOTE_MINIMAL, lineterminator="\n", *args, **kwargs
			)

	def _dict_to_list(self, rowdict):
		return super(JsonWriter, self)._dict_to_list({
			key:json.dumps(array_to_list(value))
			for key, value in rowdict.viewitems()
			})


class JsonReader(csv.DictReader, object):
	def __init__(self, *args, **kwargs):
		"""Defaults dialect=CSV_DIALECT."""
		kwargs.setdefault('dialect', CSV_DIALECT)
		super(JsonReader, self).__init__(
			quotechar = "'", quoting = csv.QUOTE_MINIMAL, *args, **kwargs
			)

		# This is a hack to strip extra quotes from the field names
		# Not proud of it, but it works.
		_ = self.fieldnames # called for side effect

		self._fieldnames = [
			fieldname.strip('"') for fieldname in self._fieldnames
			]

	def next(self):
		# type: () -> Dict[str, Any]
		attributeDict = {}  # type: Dict[Text, Any]
		for key_, raw_value_ in super(JsonReader, self).next().viewitems():
			# NOTE: Decoding UTF-8 bytes would be safer between DictReader and
			# its csv.reader as in csv32, but this is much simpler.
			if six.PY2:
				key = key_.decode('utf-8')
				raw_value = raw_value_.decode('utf-8') if raw_value_ else raw_value_
			else:
				key = key_
				raw_value = raw_value_

			try:
				value = json.loads(raw_value) if raw_value else ""

			except (ValueError, TypeError) as e:
				repr(e)
				raise Exception("failed to parse json string:{}".format(raw_value))

			match = re.search('(.*?) \((.*?)\)', key)
			if match:
				# Entry includes units so need to apply parsed units to values
				_ = units  # don't warn about `units`; it's imported for eval()
				attribute = match.group(1)
				value_units = eval(match.group(2))

				# Units do not work with empty values
				if value != '':
					# Units do not work with lists so need to convert to ndarray
					if isinstance(value, list):
						value_with_units = value_units * np.array(value)
					else:
						value_with_units = value_units * value
					# Normalize to catch any unit issues now instead of later
					value = value_with_units.normalize()

				attributeDict[attribute] = value
			else:
				attributeDict[key] = value

		# TODO(jerry): The return type differs from the base class in PY2!
		# Either suppress that mypy warning with `cast()` while we're still on
		# PY2, or go ahead and make JsonReader proxy csv.DictReader rather than
		# subclass it. Until we move to PY3 or change the base class, the type
		# checker won't know this method returns a Dict[unicode, Any].
		return cast('Dict[str, Any]', attributeDict)
