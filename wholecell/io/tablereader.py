
from __future__ import division

import os
import json
import re

from itertools import izip

import numpy as np

from . import tablewriter as tw

# TODO: tests
# TODO: load a single time point
# TODO: handle/warn/raise on inconsistent data shapes
# TODO: downsampling options

__all__ = [
	"TableReader",
	]

class TableReaderError(Exception):
	pass


class VersionError(TableReaderError):
	pass


class DoesNotExistError(TableReaderError):
	pass


class VariableWidthError(TableReaderError):
	pass


class TableReader(object):
	def __init__(self, path):

		try:
			version = open(os.path.join(path, tw.DIR_METADATA, tw.FILE_VERSION)).read()

		except IOError:
			raise VersionError("Could not open table ({}); may be wrong version".format(path))

		if version != tw.VERSION:
			raise VersionError("Expected version {} but found version {}".format(tw.VERSION, version))


		self._dirAttributes = os.path.join(path, tw.DIR_ATTRIBUTES)
		self._attributeNames = os.listdir(self._dirAttributes)

		self._dirColumns = os.path.join(path, tw.DIR_COLUMNS)
		self._columnNames = os.listdir(self._dirColumns)


	def readAttribute(self, name):
		if name not in self._attributeNames:
			raise DoesNotExistError("No such attribute: {}".format(name))

		return json.loads(
			open(os.path.join(self._dirAttributes, name)).read()
			)


	def readColumn(self, name):
		if name not in self._columnNames:
			raise DoesNotExistError("No such column: {}".format(name))

		offsets, dtype = self._loadOffsets(name)

		sizes = np.diff(offsets)

		if len(set(sizes)) > 1:
			raise VariableWidthError("Cannot load full column; data size varies")

		nEntries = sizes.size

		with open(os.path.join(self._dirColumns, name, tw.FILE_DATA)) as dataFile:
			dataFile.seek(offsets[0])

			return np.fromstring(
				dataFile.read(), dtype
				).reshape(nEntries, -1).squeeze()


	def iterColumn(self, name):
		if name not in self._columnNames:
			raise DoesNotExistError("No such column: {}".format(name))

		offsets, dtype = self._loadOffsets(name)

		sizes = np.diff(offsets)

		with open(os.path.join(self._dirColumns, name, tw.FILE_DATA)) as dataFile:
			dataFile.seek(offsets[0])

			for size in sizes:
				yield np.fromstring(
					dataFile.read(size), dtype
					)


	def readRow(self, index):
		return {
			name: self._loadEntry(name, index)
			for name in self._columnNames
			}


	def _loadEntry(self, name, index):
		offsets, dtype = self._loadOffsets(name)

		size = offsets[index+1] - offsets[index]

		with open(os.path.join(self._dirColumns, name, tw.FILE_DATA)) as dataFile:
			dataFile.seek(offsets[index])

			return np.fromstring(
				dataFile.read(size), dtype
				)


	def _loadOffsets(self, name):
		with open(os.path.join(self._dirColumns, name, tw.FILE_OFFSETS)) as offsetsFile:
			offsets = np.array([int(i.strip()) for i in offsetsFile])

		with open(os.path.join(self._dirColumns, name, tw.FILE_DATA)) as dataFile:
			rawDtype = json.loads(dataFile.read(offsets[0]))

			if isinstance(rawDtype, basestring):
				dtype = str(rawDtype)

			else:
				dtype = [ # numpy requires list-of-tuples-of-strings
					(str(n), str(t))
					for n, t in rawDtype
					]

		return offsets, dtype


	def attributeNames(self):
		return self._attributeNames


	def columnNames(self):
		return self._columnNames


	def close(self): # TODO: reimplement?
		pass
