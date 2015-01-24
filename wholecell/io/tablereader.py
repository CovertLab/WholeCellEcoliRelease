
from __future__ import division

import os
import json
import re

import numpy as np

from . import tablewriter

__all__ = [
	"TableReader",
	]


class TableReaderError(Exception):
	pass


class TableReader(object):
	def __init__(self, path):
		self._data = open(os.path.join(path, tablewriter.FILENAME_DATA))

		self._attributeDict = {}

		with open(os.path.join(path, tablewriter.FILENAME_METADATA)) as metadata:
			for line in metadata:
				name, serialized = re.split("\t", line.strip())

				value = json.loads(serialized)

				self._attributeDict[name] = value

		with open(os.path.join(path, tablewriter.FILENAME_OFFSETS)) as offsets:
			self._fieldNames = tuple(re.split("\t", offsets.readline().strip()))

			offset_values = []

			for line in offsets:
				offset_values.append(
					re.split("\t", line.strip())
					)

		self._offsets = np.array(offset_values, np.int64)

		self._nEntries = self._offsets.shape[0]

		self._startIndex = self.readAttribute(tablewriter.ATTRNAME_STARTINDEX)
		self._stepSize = self.readAttribute(tablewriter.ATTRNAME_STEPSIZE)


	def readAttribute(self, name):
		return self._attributeDict[name]


	def readColumn(self, fieldName):
		return np.array([
			self._loadData(fieldName, i) for i in xrange(self._nEntries)
			])


	def _loadData(self, fieldName, index):
		self._data.seek(
			self._offsets[index, self._fieldNames.index(fieldName)]
			)

		return np.load(self._data)


	def attributeNames(self):
		return self._attributeDict.keys()


	def fieldNames(self):
		return self._fieldNames
