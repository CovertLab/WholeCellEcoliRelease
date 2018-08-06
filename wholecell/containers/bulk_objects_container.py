
from __future__ import absolute_import
from __future__ import division

import numpy as np

class BulkObjectsContainer(object):
	"""
	A data structure for convenient, name-based operations on elements of a
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
	track the abundances (copy numbers) of molecules.  Likewise the methods
	refer to "counts" even though the data type can be floating-point.

	The number of elements and their order is inferred from the objectNames
	parameter.

	If looking at a subset of the objects repeatedly, it's best to create
	and operate on a view.

	TODO (John): Give methods more standard names.
	TODO (John): Move methods and attributes from mixedCase to under_scores.
	TODO (John): Get rid of single/group distinction, and instead check input
		types against basestring to decide what sort of output to return.
	TODO (John): Use something more generic than 'counts' to reflect the fact
		that non-integer data types are permissible.

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
		"""
		Returns the counts of all objects, or the counts associated with an
		iterable of names.

		Parameters
		----------
		names : iterable of strings
			Default is None.  If not None, the counts will be returned (in the
			provided order) for the indicated objects.

		Returns
		-------
		A vector of counts (or whatever the underlying vector represents).

		"""

		if names is None:
			return self._counts.copy()

		else:
			return self._counts[self._namesToIndexes(names)]


	def countsIs(self, values, names = None):
		"""
		Sets the counts of all objects, or the counts associated with an
		iterable of names.

		Parameters
		----------
		values : array-like
			The counts to assign to the indicated objects.
		names : iterable of strings
			Default is None.  If not None, the counts will be returned (in the
			provided order) for the indicated objects.

		"""

		if names is None:
			self._counts[:] = values

		else:
			self._counts[self._namesToIndexes(names)] = values


	def countsInc(self, values, names = None):
		"""
		Increment the counts of all objects, or the counts associated with an
		iterable of names.

		Parameters
		----------
		values : array-like
			The amounts by which to increment the counts of the indicated
			objects.
		names : iterable of strings
			Default is None.  If not None, the counts will be returned (in the
			provided order) for the indicated objects.

		"""

		values = np.asarray(values, dtype=self._counts.dtype)
		if names is None:
			self._counts[:] += values

		else:
			self._counts[self._namesToIndexes(names)] += values


	def countsDec(self, values, names = None):
		"""
		Decrement the counts of all objects, or the counts associated with an
		iterable of names.

		Parameters
		----------
		values : array-like
			The amounts by which to decrement the counts of the indicated
			objects.
		names : iterable of strings
			Default is None.  If not None, the counts will be returned (in the
			provided order) for the indicated objects.

		"""

		values = np.asarray(values, dtype=self._counts.dtype)
		if names is None:
			self._counts[:] -= values

		else:
			self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names = None):
		"""
		Returns an object that provides a permanent, ordered reference to a
		set of objects.

		Parameters
		----------
		names : iterable of strings
			Default is None.  If None, the view will simply be all objects in
			their natural order.  Otherwise the names are used to define the
			ordered elements of the view object.

		Returns
		-------
		A _BulkObjectsView instance.

		Notes
		-----
		This is the ideal way to operate on a subset of a BulkObjectsContainer
		repeatedly, as the indices only have to be gathered once.

		TODO (John): Get rid of the default behavior, as it serves no practical
			purpose or advantage over just operating on the
			BulkObjectsContainer itself.

		"""

		if names is None:
			return _BulkObjectsView(self, np.arange(self._nObjects))

		else:
			return _BulkObjectsView(self, self._namesToIndexes(names))


	def count(self, name):
		"""
		Returns the count of a single object.

		Parameters
		----------
		name : string
			The name of the indicated object.

		Returns
		-------
		The counts associated with the indicated object.

		"""
		return self._counts[self._objectIndex[name]]


	def countIs(self, value, name):
		"""
		Sets the count of a single object.

		Parameters
		----------
		value : array-like
			The count to assign to the indicated objects.
		name : string
			The name of the indicated object.

		"""
		self._counts[self._objectIndex[name]] = value


	def countInc(self, value, name):
		"""
		Increment the count of a single object.

		Parameters
		----------
		value : array-like
			The amount by which to increment the counts of the indicated
			object.
		name : string
			The name of the indicated object.

		"""
		self._counts[self._objectIndex[name]] += value


	def countDec(self, value, name):
		"""
		Decrement the count of a single object.

		Parameters
		----------
		value : array-like
			The amount by which to decrement the counts of the indicated
			object.
		name : string
			The name of the indicated object.

		"""
		self._counts[self._objectIndex[name]] -= value


	def countView(self, name):
		"""
		Returns an object that provides a permanent reference to a single
		object.

		Parameters
		----------
		name : iterable of strings
			The name of the indicated object.

		Returns
		-------
		A _BulkObjectView instance.

		Notes
		-----
		This is the ideal way to operate on an element of a
		BulkObjectsContainer repeatedly, as the index only has to be gathered
		once.

		"""
		return _BulkObjectView(self, self._objectIndex[name])

	def objectNames(self):
		"""
		Returns the names (in order) of all objects.

		Parameters
		----------
		(none)

		Returns
		-------
		A tuple of strings.

		"""
		return tuple(self._objectNames)

	def emptyLike(self):
		"""
		Returns a new BulkObjectsContainer with the same set of objects, but
		with all counts at zero.

		Notes
		-----
		This is not totally analogous to numpy.empty_like, as it fills in the
		counts with zeros (whereas numpy.empty_like's contents are arbitrary).

		"""
		names = self.objectNames()
		new_copy = BulkObjectsContainer(names)
		return new_copy

	def _namesToIndexes(self, names):
		"""
		Private method that converts an iterable of names into their
		corresponding indices into the underlying array representation.

		Parameters
		----------
		names : iterable of strings
			The names of the indicated object.

		Returns
		-------
		An arrary of indices (non-negative integers).

		"""
		return np.array([self._objectIndex[name] for name in names])


	def __eq__(self, other):
		"""
		Magic method for comparing the counts of one container to another.

		Parameters
		----------
		other : a BulkObjectContainer instance

		Returns
		-------
		True if all counts are the same, otherwise False.

		Notes
		-----
		TODO (John): This should fail if the elements of one container are
			different from another, even if their sizes are the same.

		TODO (John): If all elements are the same but in a different order,
			this method should be sensitive to that.

		TODO (John): This method really shouldn't be inspecting a private
			attribute of another object.

		"""
		return (self._counts == other._counts).all()


	def tableCreate(self, tableWriter):
		"""
		Writes the names of the objects to a 'table' file's attributes.

		Parameters
		----------
		tableWriter : a TableWriter instance

		Notes
		-----
		TODO (John): I feel that these methods pollute this class.

		"""
		tableWriter.writeAttributes(
			objectNames = self._objectNames
			)


	def tableAppend(self, tableWriter):
		"""
		Append the current counts of the objects to a 'table' file.

		Parameters
		----------
		tableWriter : a TableWriter instance

		"""
		tableWriter.append(
			counts = self._counts
			)


	def tableLoad(self, tableReader, tableIndex):
		"""
		Loads the counts of objects from a 'table' file.

		Parameters
		----------
		tableReader : a TableReader instance
		tableIndex : the row of the table from which to read the counts

		"""

		assert self._objectNames == tableReader.readAttribute("objectNames")

		self._counts = tableReader.readRow(tableIndex)["counts"]


class _BulkObjectsView(object):
	"""
	An accessor for a subset of objects in a BulkObjectsContainer.

	Parameters
	----------
	container : a BulkObjectsContainer instance
		The underlying container for the data associated with this view.
	indexes : an iterable of indices (non-negative intergers)
		The indices into the BulkObjectContainer's _counts attribute associated
		with this view.

	Notes
	-----
	TODO (John): Consider moving this class into the context of the
		BulkObjectsContainer's class definition.

	TODO (John): Consider passing the array reference rather than the container
		itself - then we don't have to access the 'private' _counts attribute.

	"""

	def __init__(self, container, indexes):
		self._container = container
		self._indexes = indexes


	def counts(self):
		"""
		Return the counts of all objects.

		Parameters
		----------
		(none)

		Returns
		-------
		A vector of counts (or whatever the underlying vector represents).

		"""
		return self._container._counts[self._indexes]


	def countsIs(self, values):
		"""
		Set the counts of all objects.

		Parameters
		----------
		values : array-like
			The counts to assign to the objects.

		"""
		self._container._counts[self._indexes] = values


	def countsInc(self, values):
		"""
		Increment the counts of all objects.

		Parameters
		----------
		values : array-like
			The amounts to add to the counts of the objects.

		"""
		values = np.asarray(values, dtype=self._container._counts.dtype)
		self._container._counts[self._indexes] += values


	def countsDec(self, values):
		"""
		Decrement the counts of all objects.

		Parameters
		----------
		values : array-like
			The amounts to subtract from the counts of the objects.

		"""
		values = np.asarray(values, dtype=self._container._counts.dtype)
		self._container._counts[self._indexes] -= values


class _BulkObjectView(object):
	"""
	An accessor for a single object in a BulkObjectsContainer.

	Parameters
	----------
	container : a BulkObjectsContainer instance
		The underlying container for the data associated with this view.
	index : a non-negative integer
		The index into the BulkObjectContainer's _counts attribute associated
		with this view.

	Notes
	-----
	TODO (John): Consider moving this class into the context of the
		BulkObjectsContainer's class definition.

	TODO (John): Consider passing the array reference rather than the container
		itself - then we don't have to access the 'private' _counts attribute.

	"""

	def __init__(self, container, index):
		self._container = container
		self._index = index


	def count(self):
		"""
		Return the count of the object.

		Parameters
		----------
		(none)

		Returns
		-------
		The counts associated with the indicated object.

		"""
		return self._container._counts[self._index]


	def countIs(self, values):
		"""
		Set the count of the objects.

		Parameters
		----------
		value : array-like
			The count to assign to the object.

		"""
		self._container._counts[self._index] = values


	def countInc(self, values):
		"""
		Increment the counts of the object.

		Parameters
		----------
		value : array-like
			The amount to add to the count of the object.

		"""
		self._container._counts[self._index] += values


	def countDec(self, values):
		"""
		Decrement the counts of the object.

		Parameters
		----------
		value : array-like
			The amount to subtract from the count of the object.

		"""
		self._container._counts[self._index] -= values
