"""
Subclasses of DictWriter and DictReader that parse plaintext as JSON strings,
allowing for basic type parsing and fields that are dictionaries or lists.
"""

import csv
import json
import re
import numpy as np

from wholecell.utils import units

def array_to_list(value):
	if isinstance(value, np.ndarray):
		value = value.tolist()

	return value


class JsonWriter(csv.DictWriter):
	def __init__(self, *args, **kwargs):
		csv.DictWriter.__init__(
			self, quotechar = "'", quoting = csv.QUOTE_MINIMAL, lineterminator="\n", *args, **kwargs
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
		attributeDict = {}
		for key, raw_value in csv.DictReader.next(self).viewitems():
			try:
				value = json.loads(raw_value) if raw_value else ""

			except (ValueError, TypeError) as e:
				repr(e)
				raise Exception("failed to parse json string:{}".format(raw_value))

			try:
				attribute = re.search('(.*?) \(', key).group(1)
				value_units =  eval(re.search('\((.*?)\)',key).group(1))
				attributeDict[attribute] = value * value_units
			except AttributeError:
				attributeDict[key] = value
		return attributeDict

		# return {
		# 	key:json.loads(value) if value else "" # catch for empty field
		# 	for key, value in csv.DictReader.next(self).viewitems()
		# 	}
