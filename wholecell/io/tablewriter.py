
from __future__ import division

import os
import json

import numpy as np

__all__ = [
	"TableWriter",
	# "TableWriterError",
	# "FilesClosedError",
	# "MissingFieldError",
	# "UnrecognizedFieldError",
	# "AttributeAlreadyExistsError",
	# "AttributeTypeError"
	]

FILENAME_DATA = "data"
FILENAME_METADATA = "metadata"
FILENAME_OFFSETS = "offsets"

ATTRNAME_STARTINDEX = "startIndex"
ATTRNAME_STEPSIZE = "stepSize"

class TableWriterError(Exception):
	pass

class FilesClosedError(TableWriterError):
	pass

class MissingFieldError(TableWriterError):
	pass

class UnrecognizedFieldError(TableWriterError):
	pass

class AttributeAlreadyExistsError(TableWriterError):
	pass

class AttributeTypeError(TableWriterError):
	pass

class TableWriter(object):
	def __init__(self, path, stepSize = 1, startIndex = 0):

		os.mkdir(path)

		self._data = open(os.path.join(path, FILENAME_DATA), "w")
		self._metadata = open(os.path.join(path, FILENAME_METADATA), "w")
		self._offsets = open(os.path.join(path, FILENAME_OFFSETS), "w")

		self._stepSize = stepSize
		self._index = startIndex

		self._fieldNames = None
		self._attributeNames = []

		self._closed = False

		self.writeAttributes(**{
			ATTRNAME_STARTINDEX:self._index,
			ATTRNAME_STEPSIZE:self._stepSize
			})


	def append(self, **fieldsAndValues):
		if self._closed:
			raise FilesClosedError()

		if self._fieldNames is None:
			self._fieldNames = tuple(sorted(fieldsAndValues.keys()))

			self._offsets.write("\t".join(self._fieldNames) + "\n")

		else:
			missingFields = set(self._fieldNames) - fieldsAndValues.viewkeys()
			unrecognizedFields = fieldsAndValues.viewkeys() - set(self._fieldNames)

			if missingFields:
				raise MissingFieldError(
					"Missing fields: {}".format(missingFields)
					)

			if unrecognizedFields:
				raise UnrecognizedFieldError(
					"Unrecognized fields: {}".format(unrecognizedFields)
					)

		offsets = []

		for field in self._fieldNames:
			offsets.append(self._data.tell())

			np.save(self._data, fieldsAndValues[field])

		self._offsets.write(
			"\t".join("{}".format(offset) for offset in offsets) + "\n"
			)

		self._index += self._stepSize


	def writeAttributes(self, **namesAndValues):
		if self._closed:
			raise FilesClosedError()

		for name, value in namesAndValues.viewitems():
			if name in self._attributeNames:
				raise AttributeAlreadyExistsError(
					"An attribute named {} already exists.".format(name)
					)

			try:
				serialized = json.dumps(value)

			except TypeError:
				raise AttributeTypeError(
					"Attribute {} value ({}) was not JSON serializable.".format(
						name, value # TODO: repr for value instead of str
						)
					)

			self._metadata.write(
				"{}\t{}\n".format(name, serialized)
				)

			self._attributeNames.append(name)


	def close(self):
		if not self._closed:
			self._data.close()
			self._metadata.close()
			self._offsets.close()

			self._closed = True
