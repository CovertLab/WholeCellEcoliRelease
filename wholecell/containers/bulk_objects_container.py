
from __future__ import absolute_import
from __future__ import division

import numpy as np

class BulkObjectsContainer(object):
	"""
	A data structure for convenient, name-based operations on elements of a
	NumPy vector.

	The BulkObjectsContainer provides an interface to the elements of a NumPy
	vector by name.  Oftentimes in a multi-domain, multi-paradigm simulation,
	we wish to store all model variables in one place, yet will only interact
	with a subset of variables in each submodel.  This class makes it easy to
	access and modify those variable values.  For example:

	>>> names = ['A', 'B', 'C', 'D'] # names of some model variables
	>>> container = BulkObjectsContainer(names)
	>>> container.countIs('A', 10) # set A to 10
	>>> container.countsInc(5) # increase the counts of everything by 5
	>>> container.counts(['D', 'C']) # access D and C, in that order

	You can also create a 'view' to interact with the data in a fixed order:

	>>> view = container.countsView(['D', 'C'])
	>>> view.countsIs([10, 15]) # set D to 10 and C to 15
	>>> view.countsDec([3, 7]) # decrease D by 3 and C by 7
	>>> view.counts() # get the counts of D and C

	If used more than once, creating a view is more efficient than operating on
	the container directly.

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
	The default data type is integers because the original use case was to
	track the abundances (copy numbers) of molecules.  Likewise the methods
	refer to "counts" even though the data type can be floating-point.

	Floating-point use cases:
	- concentrations of molecules
	- mean or standard deviations of counts

	The number of elements and their order is inferred from the objectNames
	parameter.

	TODO (John): Give methods more standard names (e.g. set, get, add).
	TODO (John): Move methods and attributes from mixedCase to under_scores.
	TODO (John): Get rid of single/group distinction in methods/views, and
		instead check input types against basestring to decide what sort of
		output to return.
	TODO (John): Use something more generic than 'counts' to reflect the fact
		that non-integer data types are permissible.

	"""

	def __init__(self, objectNames, dtype = np.int64):
		# Copy the object names into a tuple to ensure they are ordered and
		# immutable
		self._objectNames = tuple(objectNames)
		self._nObjects = len(self._objectNames)

		self._dtype = dtype

		# Store the indices for each element in a dictionary for faster
		# (on average, O(1)) look-up (list.index is slow, O(n))
		self._objectIndex = {
			objectName:index
			for index, objectName in enumerate(self._objectNames)
			}

		self._counts = np.zeros(self._nObjects, self._dtype)


	def counts(self, names = None):
		"""
		Get the counts of objects.

		Parameters
		----------
		names : iterable of strings
			The names of the objects.  If None (default), then all objects are
			used in their original ordering.

		Returns
		-------
		A vector (1D numpy.ndarray) of counts.

		"""

		if names is None:
			return self._counts.copy()

		else:
			return self._counts[self._namesToIndexes(names)]


	def countsIs(self, values, names = None):
		"""
		Set the counts of objects.

		Parameters
		----------
		values : array-like
			The assigned counts.
		names : iterable of strings
			The names of the objects.  If None (default), then all objects are
			used in their original ordering.

		"""

		if names is None:
			self._counts[:] = values

		else:
			self._counts[self._namesToIndexes(names)] = values


	def countsInc(self, values, names = None):
		"""
		Increment the counts of objects.

		Parameters
		----------
		values : array-like
			The added counts.
		names : iterable of strings
			The names of the objects.  If None (default), then all objects are
			used in their original ordering.

		"""

		values = np.asarray(values, dtype=self._dtype)
		if names is None:
			self._counts[:] += values

		else:
			self._counts[self._namesToIndexes(names)] += values


	def countsDec(self, values, names = None):
		"""
		Decrement the counts of objects.

		Parameters
		----------
		values : array-like
			The subtracted counts.
		names : iterable of strings
			The names of the objects.  If None (default), then all objects are
			used in their original ordering.

		"""

		values = np.asarray(values, dtype=self._dtype)
		if names is None:
			self._counts[:] -= values

		else:
			self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names = None):
		"""
		Create a view on a set of objects.

		Parameters
		----------
		names : iterable of strings
			The names of the objects.  If None (default), then all objects are
			used in their original ordering.

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
		Get the count of an object.

		Parameters
		----------
		name : object name

		Returns
		-------
		A scalar.

		"""
		return self._counts[self._objectIndex[name]]


	def countIs(self, value, name):
		"""
		Set the count of an object.

		Parameters
		----------
		value : scalar
			The assigned count.
		name : object name

		"""
		self._counts[self._objectIndex[name]] = value


	def countInc(self, value, name):
		"""
		Increment the count of an object.

		Parameters
		----------
		value : scalar
			The added count.
		name : object name

		"""
		self._counts[self._objectIndex[name]] += value


	def countDec(self, value, name):
		"""
		Decrement the count of an object.

		Parameters
		----------
		value : scalar
			The subtracted count.
		name : object name

		"""
		self._counts[self._objectIndex[name]] -= value


	def objectNames(self):
		"""
		The names (in order) of all objects.

		Parameters
		----------
		(none)

		Returns
		-------
		A tuple of strings.

		"""
		return self._objectNames

	def emptyLike(self):
		"""
		Create a new BulkObjectsContainer with the same set of objects, but
		with all counts at zero.

		Notes
		-----
		This is not totally analogous to numpy.empty_like, as it fills in the
		counts with zeros (whereas numpy.empty_like's contents are arbitrary).

		"""
		names = self.objectNames()
		new_copy = BulkObjectsContainer(names, dtype = self._dtype)
		return new_copy

	def _namesToIndexes(self, names):
		"""
		Convert an iterable of names into their corresponding indices into the
		underlying array representation.

		Parameters
		----------
		names : iterable of strings
			The names of the objects.

		Returns
		-------
		An arrary of indices (non-negative integers).

		TODO (John): Handle the case that names is None, returning
			np.arange(self._nObjects).  This would simplify many other methods.
			We can't use slice(None) because we need the data to be copied on
			retrieval.  Otherwise modifying an output array (e.g. an array
			return by .counts()) would modify the underlying data.

		"""
		return np.array([self._objectIndex[name] for name in names])


	def __eq__(self, other):
		"""
		Return whether the counts of one BulkObjectsContainer instance are all
		equal to another.

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
		Write the names of the objects to a 'table' file's attributes.

		Parameters
		----------
		tableWriter : a TableWriter instance

		Notes
		-----
		TODO (John): I feel that these methods pollute this class.  This
			interface is simple enough that we could remove it entirely, or
			we could add functions to this module that handle the interface
			between TableWriter and BulkObjectsContainer without polluting
			either class with information about the other.

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
		Load the counts of objects from a 'table' file.

		Parameters
		----------
		tableReader : a TableReader instance
		tableIndex : non-negative integer
			Specify which row of the table to row.

		TODO (John): If this was a class method, it could instantiate with the
			correct object names instead of asserting that the names are
			correct.

		"""

		assert self._objectNames == tuple(tableReader.readAttribute("objectNames"))

		self._counts = tableReader.readRow(tableIndex)["counts"]


class _BulkObjectsView(object):
	"""
	A view onto selected objects in a BulkObjectsContainer.

	Parameters
	----------
	container : a BulkObjectsContainer instance
	indexes : an iterable of indices (non-negative integers)
		Specify the indices into the BulkObjectContainer's _counts attribute
		associated with this view's objects.

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
		Get the counts of all objects.

		Parameters
		----------
		(none)

		Returns
		-------
		A vector (1D numpy.ndarray) of counts.

		"""
		return self._container._counts[self._indexes]


	def countsIs(self, values):
		"""
		Set the counts of all objects.

		Parameters
		----------
		values : array-like
			The assigned counts.

		"""
		self._container._counts[self._indexes] = values


	def countsInc(self, values):
		"""
		Increment the counts of all objects.

		Parameters
		----------
		values : array-like
			The added counts.

		"""
		values = np.asarray(values, dtype=self._container._dtype)
		self._container._counts[self._indexes] += values


	def countsDec(self, values):
		"""
		Decrement the counts of all objects.

		Parameters
		----------
		values : array-like
			The subtracted counts.

		"""
		values = np.asarray(values, dtype=self._container._dtype)
		self._container._counts[self._indexes] -= values
