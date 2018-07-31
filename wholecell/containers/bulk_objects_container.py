'''

bulk_objects_container.py

'''

from __future__ import absolute_import
from __future__ import division

import numpy as np


class BulkObjectsContainer(object):
	"""
	An data structure for convenient, name-based operations on elements of a
	NumPy vector.

	Parameters
	----------
	objectNames : iterable of strings
		The names that will be used to refer to the elements of the
		underlying vector.
	dtype : a valid NumPy datatype identifier (default: np.int64)
		The data type of the underlying array.

	See also
	--------
	wholecell.containers.unique_objects_container.UniqueObjectsContainer

	Notes
	-----
	The dafault data type is integers because the original use case was to
	track the abundances (copy numbers) of molecules.

	The number of elements and their order is inferred from the objectNames
	parameter.

	TODO (John): Give methods more standard names.
	TODO (John): Move methods and attributes from mixedCase to under_scores.

	"""

	def __init__(self, objectNames, dtype = np.int64):
		# Cast objectNames to list to make sure it is ordered, and to prevent
		# any side effects from external code modifying the parameter
		self._objectNames = list(objectNames)

		self._nObjects = len(self._objectNames)

		# Store the indices for each element in a dictionary for faster
		# look-up (list.index is slow)
		self._objectIndex = {
			objectName:index
			for index, objectName in enumerate(self._objectNames)
			}

		self._counts = np.zeros(len(self._objectNames), dtype)


	def counts(self, names = None):
		if names is None:
			return self._counts.copy()

		else:
			return self._counts[self._namesToIndexes(names)]


	def countsIs(self, values, names = None):
		if names is None:
			self._counts[:] = values

		else:
			self._counts[self._namesToIndexes(names)] = values


	def countsInc(self, values, names = None):
		values = np.asarray(values, dtype=self._counts.dtype)
		if names is None:
			self._counts[:] += values

		else:
			self._counts[self._namesToIndexes(names)] += values


	def countsDec(self, values, names = None): # TODO: raise exception if > max?
		values = np.asarray(values, dtype=self._counts.dtype)
		if names is None:
			self._counts[:] -= values

		else:
			self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names = None):
		if names is None:
			return _BulkObjectsView(self, np.arange(self._nObjects))

		else:
			return _BulkObjectsView(self, self._namesToIndexes(names))


	def count(self, name):
		return self._counts[self._objectIndex[name]]


	def countIs(self, value, name):
		self._counts[self._objectIndex[name]] = value


	def countInc(self, value, name):
		self._counts[self._objectIndex[name]] += value


	def countDec(self, value, name):
		self._counts[self._objectIndex[name]] -= value


	def countView(self, name):
		return _BulkObjectView(self, self._objectIndex[name])

	def objectNames(self):
		return tuple(self._objectNames)

	def emptyLike(self):
		names = self.objectNames()
		new_copy = BulkObjectsContainer(names)
		return new_copy

	def _namesToIndexes(self, names):
		return np.array([self._objectIndex[name] for name in names])


	def __eq__(self, other):
		return (self._counts == other._counts).all()


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			objectNames = self._objectNames
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			counts = self._counts
			)


	def tableLoad(self, tableReader, tableIndex):
		assert self._objectNames == tableReader.readAttribute("objectNames")

		self._counts = tableReader.readRow(tableIndex)["counts"]


class _BulkObjectsView(object):
	'''
	_BulkObjectsView

	An accessor for a subset of objects in a BulkObjectsContainer.
	'''

	def __init__(self, container, indexes):
		self._container = container
		self._indexes = indexes


	def counts(self):
		return self._container._counts[self._indexes]


	def countsIs(self, values):
		self._container._counts[self._indexes] = values


	def countsInc(self, values):
		values = np.asarray(values, dtype=self._container._counts.dtype)
		self._container._counts[self._indexes] += values


	def countsDec(self, values):
		values = np.asarray(values, dtype=self._container._counts.dtype)
		self._container._counts[self._indexes] -= values


class _BulkObjectView(object):
	'''
	_BulkObjectView

	An accessor for a single object in a BulkObjectsContainer.
	'''

	def __init__(self, container, index):
		self._container = container
		self._index = index


	def count(self):
		return self._container._counts[self._index]


	def countIs(self, values):
		self._container._counts[self._index] = values


	def countInc(self, values):
		self._container._counts[self._index] += values


	def countDec(self, values):
		self._container._counts[self._index] -= values
