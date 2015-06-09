
from __future__ import division

import os
import json

import numpy as np

# TODO: tests

__all__ = [
	"TableWriter",
	# "TableWriterError",
	# "FilesClosedError",
	# "MissingFieldError",
	# "UnrecognizedFieldError",
	# "AttributeAlreadyExistsError",
	# "AttributeTypeError"
	]

VERSION = "2" # should update this any time there is a spec-breaking change

DIR_METADATA = "metadata"
DIR_ATTRIBUTES = "attributes"
DIR_COLUMNS = "columns"
FILE_VERSION = "version"
FILE_DATA = "data"
FILE_OFFSETS = "offsets"

def _silent_makedirs(path):
	try:
		os.makedirs(path)

	except OSError: # thrown if folders exist
		pass

	return path

class TableWriterError(Exception):
	pass

class MissingFieldError(TableWriterError):
	pass

class UnrecognizedFieldError(TableWriterError):
	pass

class AttributeAlreadyExistsError(TableWriterError):
	pass

class AttributeTypeError(TableWriterError):
	pass


class _Column(object):
	def __init__(self, path):
		_silent_makedirs(path)

		self._data = open(os.path.join(path, FILE_DATA), "w")
		self._offsets = open(os.path.join(path, FILE_OFFSETS), "w")

		self._dtype = None


	def append(self, value):
		value = np.asarray(value, self._dtype)

		if self._dtype is None:
			self._dtype = value.dtype

			descr = self._dtype.descr
			if len(descr) == 1 and descr[0][0] == "":
				descr = descr[0][1]

			self._data.write(json.dumps(descr) + "\n")
			self._offsets.write(str(self._data.tell()) + "\n")

		self._data.write(value.tobytes())
		self._offsets.write(str(self._data.tell()) + "\n")


	def close(self):
		self._data.close()
		self._offsets.close()


	def __del__(self):
		self.close()


class TableWriter(object):
	def __init__(self, path):

		dirMetadata = _silent_makedirs(os.path.join(path, DIR_METADATA))

		open(os.path.join(dirMetadata, FILE_VERSION), "w").write(VERSION)

		self._dirAttributes = _silent_makedirs(os.path.join(path, DIR_ATTRIBUTES))
		self._attributeNames = []

		self._dirColumns = _silent_makedirs(os.path.join(path, DIR_COLUMNS))
		self._columns = None


	def append(self, **namesAndValues):
		if self._columns is None:
			self._columns = {
				name:_Column(os.path.join(self._dirColumns, name))
				for name in namesAndValues.viewkeys()
				}

		else:
			missingFields = self._columns.viewkeys() - namesAndValues.viewkeys()
			unrecognizedFields = namesAndValues.viewkeys() - self._columns.viewkeys()

			if missingFields:
				raise MissingFieldError(
					"Missing fields: {}".format(", ".join(missingFields))
					)

			if unrecognizedFields:
				raise UnrecognizedFieldError(
					"Unrecognized fields: {}".format(", ".join(unrecognizedFields))
					)

		for name, value in namesAndValues.viewitems():
			self._columns[name].append(value)


	def writeAttributes(self, **namesAndValues):
		for name, value in namesAndValues.viewitems():
			if name in self._attributeNames:
				raise AttributeAlreadyExistsError(
					"An attribute named '{}' already exists.".format(name)
					)

			try:
				if isinstance(value, np.ndarray):
					print "Warning - converting '{}' attribute from ndarray to list for JSON serialization.".format(name)

					value = value.tolist()

				serialized = json.dumps(value)

			except TypeError:

				raise AttributeTypeError(
					"Attribute '{}' value ({}) was not JSON serializable.".format(
						name, value # TODO: repr for value instead of str
						)
					)

			open(os.path.join(self._dirAttributes, name), "w").write(serialized)

			self._attributeNames.append(name)


	def close(self):
		if self._columns is not None:
			for column in self._columns.viewvalues():
				column.close()


	def __del__(self):
		self.close()
