'''

bulk_objects_container.py

'''

import numpy as np


class BulkObjectsContainer(object):
	'''
	BulkObjectsContainer

	A wrapper around a NumPy array that tracks the "bulk" counts of objects.
	Bulk objects are those which have no identity beyond their names.
	'''

	def __init__(self, objectNames):
		self._objectNames = objectNames

		self._objectIndex = {objectName:index for index, objectName in enumerate(objectNames)}

		self._counts = np.zeros(len(objectNames), np.uint64)


	def counts(self, names = None):
		if names == None:
			return self._counts.copy()

		else:
			return self._counts[self._namesToIndexes(names)]


	def countsIs(self, names, values):
		self._counts[self._namesToIndexes(names)] = values


	def countsInc(self, names, values): 
		self._counts[self._namesToIndexes(names)] += values


	def countsDec(self, names, values): # TODO: raise exception if > max?
		self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names):
		raise _BulkObjectsView(self, self._namesToIndexes(names))


	def count(self, name):
		return self._counts[self._objectIndex[name]]


	def countInc(self, name, value):
		self._counts[self._objectIndex[name]] += value


	def countDec(self, name, value):
		self._counts[self._objectIndex[name]] -= value


	def countView(self, name):
		raise _BulkObjectView(self, self._objectIndex[name])


	def _namesToIndexes(self, names):
		return np.array([self._objectIndex[name] for name in names])

	# TODO: mass calculation
	# TODO: saving


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
		self._container._counts[self._indexes] += values


	def countsDec(self, values):
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
