"""
Subclasses of DictWriter and DictReader that parse plaintext as JSON strings,
allowing for basic type parsing and fields that are dictionaries or lists.
"""

import csv
import json

import numpy as np

def array_to_list(value):
	if isinstance(value, np.ndarray):
		value = value.tolist()

	return value


class JsonWriter(csv.DictWriter):
	def __init__(self, *args, **kwargs):
		csv.DictWriter.__init__(
			self, quotechar = "'", quoting = csv.QUOTE_MINIMAL, *args, **kwargs
			)

	def _dict_to_list(self, rowdict):
		return csv.DictWriter._dict_to_list(self, {
			key:json.dumps(array_to_list(value))
			for key, value in rowdict.viewitems()
			})


class JsonReader(csv.DictReader):
	def __init__(self, *args, **kwargs):
		csv.DictReader.__init__(
			self, quotechar = "'", quoting = csv.QUOTE_MINIMAL, *args, **kwargs
			)

		# This is a hack to strip extra quotes from the field names
		# Not proud of it, but it works.

		self.fieldnames # called for side effect

		self._fieldnames = [
			fieldname.strip('"') for fieldname in self._fieldnames
			]

	def next(self):
		return {
			key:json.loads(value) if value else "" # catch for empty field
			for key, value in csv.DictReader.next(self).viewitems()
			}
