"""
A UniqueObjectsContainer tracks the attributes of unique objects, which are
typically molecules.

The UniqueObjectsContainer uses _UniqueObject instances to present O-O proxies
to specific molecules.
"""

from __future__ import absolute_import, division, print_function

from copy import deepcopy
from itertools import izip
from functools import partial

import numpy as np
import zlib

import wholecell.utils.linear_programming as lp

# TODO: object transfer between UniqueObjectsContainer instances
# TODO: unique id for each object based on
#	hash(time, collectionIndex, arrayIndex, containerName) at creation time

_MAX_ID_SIZE = 40 # max length of the unique id assigned to objects

class UniqueObjectsContainerException(Exception):
	pass


def decomp(specifications, compressed_collections, global_ref_count):
	"""Decompress the arguments into a UniqueObjectsContainer. "decomp" is
	intentionally short and awkward for intended use limited to pickling. It
	calls the constructor to set up indexes and caches, unlike `__setstate__`.

	CAUTION: Future edits are expected to maintain backward compatibility with
	stored pickled arguments (e.g. via optional args) or explicitly detach
	(e.g. `__reduce__` to a new unpickling function.)

	The current stored format is sensitive to all the internal representation
	details including inactive entries, bookkeeping fields, and index fields
	which could all be reconstructed here.

	Args:
		specifications (dict): dtype specs as passed to UniqueObjectsContainer()
		compressed_collections (list[bytes]): zlib-compressed bytes of the collections ndarrays
		global_ref_count (int): the size of the global references ndarray

	Returns:
		A filled-in UniqueObjectsContainer.
	"""
	container = UniqueObjectsContainer(specifications)
	_collections = container._collections

	# Decompress the _collections arrays.
	for index, data in enumerate(compressed_collections):
		_collections[index] = decompress_ndarray(data, _collections[index].dtype)

	# Reconstruct the _globalReference array.
	grefs = np.zeros(global_ref_count, container._globalReference.dtype)
	container._globalReference = grefs
	for collectionIndex, collection in enumerate(_collections):
		for objectIndex, entry in enumerate(collection):
			if entry["_entryState"] == container._entryActive:
				global_ref = grefs[entry["_globalIndex"]]
				global_ref["_entryState"] = container._entryActive
				global_ref["_collectionIndex"] = collectionIndex
				global_ref["_objectIndex"] = objectIndex

	return container


def compress_ndarray(array):
	return zlib.compress(array.tobytes(), 9)

def decompress_ndarray(compressed_bytes, dtype):
	return np.frombuffer(zlib.decompress(compressed_bytes), dtype)


def make_dtype_spec(dict_of_dtype_specs):
	"""Construct a structured dtype spec from a dict of {name: dtype} specs.
	Sort it to ensure a consistent storage format for the structured values.
	"""
	return sorted([
		(attrName, attrType)
		for attrName, attrType in dict_of_dtype_specs.iteritems()
		])


class UniqueObjectsContainer(object):
	"""
	Essentially a database of unique molecules and other unique objects kept in
	a dict of structured arrays (that is, a dict of ndarrays of structs). Each
	structured array (DB table) names a collection of similar objects, each
	array entry (DB row) holds the state for a unique object instance, and the
	structured array fields (DB columns) hold its attributes.
	Used for unique molecules state and partitions.

	Parameters:
		specifications (Dict[str, Dict[str, str]]): Maps the unique molecule
			names (collection names) to the {attribute_name: dtype} molecule
			attributes (structured array fields). The dtype declarations are
			strings describing NumPy scalar types.

			Example: {
				'DNA polymerase': {'bound': 'bool', 'location': 'int32'},
				'RNA polymerase': {'bound': 'bool', 'location': 'int32', 'decay': 'float32'},
				}

	You can store a UniqueObjectsContainer via TableWriter or (for a single
	snapshot) more compactly via pickling.
	"""

	# State descriptions
	_entryInactive = 0 # an available entry; 0 works for np.zeros() allocation
	_entryActive = 1 # an entry that is in use
	# TODO(jerry): Store the _entryState in a parallel array? Or keep a count
	#   of active entries and always compact them? Do int8 _entryState fields
	#   make in-memory performance faster than int64 (better cache locality) or
	#   slower (maybe misaligned following fields)? What about in TableWriter?
	#   When pickling, the data gets compressed so the impact is smaller.
	# TODO(jerry): Use narrower index fields?

	_defaultSpecification = {  # bookkeeping fields to add to every struct type
		"_entryState":np.int8, # see state descriptions above
		"_globalIndex":np.int64, # index in the _globalReference array (collection)
		# "_uniqueId":"{}str".format(_MAX_ID_SIZE), # unique ID assigned to each object
		}

	_globalReferenceDtype = {
		"_entryState":np.int8, # see state descriptions above
		"_collectionIndex":np.int64,
		"_objectIndex":np.int64,
		}

	_fractionExtendEntries = 0.1 # fractional rate to increase number of entries in the structured array (collection)

	_queryOperations = {
		">":np.greater,
		">=":np.greater_equal,
		"<":np.less,
		"<=":np.less_equal,
		"==":np.equal,
		"!=":np.not_equal,
		"in":np.lib.arraysetops.in1d,
		"not in":partial(np.lib.arraysetops.in1d, invert = True)
		}

	def __init__(self, specifications):
		self._collections = [] # ordered list of numpy structured arrays
		self._nameToIndexMapping = {} # collectionName:index of associated structured array

		self._specifications = deepcopy(specifications) # collectionName:{attributeName:type}

		self._names = tuple(sorted(self._specifications.keys())) # sorted collection names

		defaultSpecKeys = self._defaultSpecification.viewkeys()

		# Add the attributes used internally
		for name, specification in self._specifications.viewitems():
			# Make sure there is no overlap in attribute names
			invalidAttributeNames = (specification.viewkeys() & defaultSpecKeys)
			if invalidAttributeNames:
				raise UniqueObjectsContainerException(
					"Invalid attribute names in specification for {}: {}".format(
						name,
						', '.join(invalidAttributeNames)
						)
					)

			specification.update(self._defaultSpecification)

		for collectionIndex, collectionName in enumerate(self._names):
			specification = self._specifications[collectionName]

			# Create the structured array collection with 1 inactive entry.
			newArray = np.zeros(1, dtype = make_dtype_spec(specification))

			# Create references to collections
			self._collections.append(newArray)
			self._nameToIndexMapping[collectionName] = collectionIndex

		# Create an array which handles global references to objects in all
		# _collections. This is an index across all the unique object names
		# (molecule types) in this container.
		self._globalReference = np.zeros(
			1, dtype = make_dtype_spec(self._globalReferenceDtype))


	def __reduce__(self):
		"""Reduce the container to its defining state for pickling.
		Compress the state for transmission efficiency.
		Return a callable object and its args.
		"""
		# TODO(jerry): Compress the specs?
		# TODO(jerry): Squeeze out inactive entries and _entryState fields, but
		#   indicate which slots are inactive or else make __eq__ ignore
		#   differences in their positions and impact on indexes. Worth the work?
		specs = self._copy_specs()
		compressed_collections = [compress_ndarray(col) for col in self._collections]
		return decomp, (specs, compressed_collections, self._globalReference.size)


	def _growArray(self, array, nObjects):
		"""Grow the array if needed to make room for nObjects, then return the
		new array and the indexes of nObjects inactive entries.
		"""
		freeIndexes = np.where(
			array["_entryState"] == self._entryInactive
			)[0]

		nFreeIndexes = freeIndexes.size

		newArray = array
		if nFreeIndexes < nObjects:
			oldSize = array.size
			nNewEntries = max(
				np.int64(oldSize * self._fractionExtendEntries),
				nObjects - nFreeIndexes
				)

			newArray = np.append(
				array,
				np.zeros(nNewEntries, dtype = array.dtype)
				)

			freeIndexes = np.concatenate((
				freeIndexes,
				oldSize + np.arange(nNewEntries)
				))

		return newArray, freeIndexes[:nObjects]

	def _getFreeIndexes(self, collectionIndex, nObjects):
		"""Return indexes of nObjects inactive entries in the specified
		collection and in _globalReference after extending the arrays if needed.
		"""
		self._collections[collectionIndex], freeIndexes = self._growArray(
			self._collections[collectionIndex], nObjects)
		self._globalReference, freeGlobalIndexes = self._growArray(
			self._globalReference, nObjects)
		return freeIndexes, freeGlobalIndexes


	def objectsNew(self, collectionName, nObjects, **attributes):
		"""Add nObjects new objects/molecules of the named type, all with the
		given attributes. Returns a _UniqueObjectSet proxy for the new entries.
		"""
		collectionIndex = self._nameToIndexMapping[collectionName]
		objectIndexes, globalIndexes = self._getFreeIndexes(collectionIndex, nObjects)

		collection = self._collections[collectionIndex]

		# TODO: restore unique object IDs
		# TODO(jerry): Would it be faster to copy one new entry to all rows
		# then set the _globalIndex columns?

		collection["_entryState"][objectIndexes] = self._entryActive
		collection["_globalIndex"][objectIndexes] = globalIndexes
		# collection["_uniqueId"][objectIndexes] = uniqueObjectIds

		for attrName, attrValue in attributes.viewitems():
			collection[attrName][objectIndexes] = attrValue

		self._globalReference["_entryState"][globalIndexes] = self._entryActive
		self._globalReference["_collectionIndex"][globalIndexes] = collectionIndex
		self._globalReference["_objectIndex"][globalIndexes] = objectIndexes

		return _UniqueObjectSet(self, globalIndexes)


	def objectNew(self, collectionName, **attributes):
		"""Add a new object/molecule of the named type with the given
		attributes. Returns a _UniqueObject proxy for the new entry.
		"""
		(molecule,) = self.objectsNew(collectionName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def objectsDel(self, objects):
		"""Delete unique objects given a _UniqueObjectSet or other iterable of
		_UniqueObject proxies.
		"""
		for obj in objects:
			self.objectDel(obj)


	def objectDel(self, obj):
		"""Delete a unique object given its _UniqueObject proxy."""
		# TODO: profile this in processes that delete lots of objects
		# TODO(jerry): This should just call a public method on _UniqueObject
		# that modifies the proxy and calls a private method to modify the
		# container, which could be a different container.
		assert obj._container == self

		collection = self._collections[obj._collectionIndex]
		collection[obj._objectIndex].fill(0)

		self._globalReference[obj._globalIndex].fill(0)

		obj._globalIndex = -1
		obj._collectionIndex = -1
		obj._objectIndex = -1


	def objects(self, **operations):
		"""Return a _UniqueObjectSet proxy for all objects (molecules) that
		satisfy an optional attribute query. Querying every object is generally
		not what you want to do. The queried attributes must be in all the
		collections (all molecule types).
		"""
		if operations:
			results = []

			for collectionIndex in xrange(len(self._collections)):
				results.append(self._queryObjects(collectionIndex, **operations))

			return _UniqueObjectSet(self, np.concatenate([
				self._collections[collectionIndex]["_globalIndex"][result]
				for collectionIndex, result in enumerate(results)
				]))

		else:
			return _UniqueObjectSet(self,
				np.where(self._globalReference["_entryState"] == self._entryActive)[0]
				)


	def objectsInCollection(self, collectionName, **operations):
		"""Return a _UniqueObjectSet proxy for all objects (molecules) belonging
		to a named collection that satisfy an optional attribute query.
		"""
		# TODO(jerry): Special case the empty query?
		collectionIndex = self._nameToIndexMapping[collectionName]

		result = self._queryObjects(collectionIndex, **operations)

		return _UniqueObjectSet(self,
			self._collections[collectionIndex]["_globalIndex"][result]
			)


	def objectsInCollections(self, collectionNames, **operations):
		"""Return a _UniqueObjectSet proxy for all objects (molecules)
		belonging to the given collection names that satisfy an optional
		attribute query. The queried attributes must be in all the named
		collections.
		"""
		collectionIndexes = [self._nameToIndexMapping[collectionName] for collectionName in collectionNames]
		results = []

		for collectionIndex in collectionIndexes:
			results.append(self._queryObjects(collectionIndex, **operations))

		return _UniqueObjectSet(self, np.concatenate([
			self._collections[collectionIndex]["_globalIndex"][result]
			for collectionIndex, result in izip(collectionIndexes, results)
			]))


	def _queryObjects(self, collectionIndex, **operations):
		"""For the structured array at the given index, perform the given
		comparison operations and return a boolean array that's True where
		entries passed all the comparisons. `operations` is a dict
		`{attribute_name: (comparison_operator_str, query_value)}`.

		The comparison_operator_str is one of '>', '>=', '<', '<=', '==', '!=',
		'in', 'not in'. `operations` can perform at most one comparison per
		attribute. The matched sets are intersected together.
		"""
		operations["_entryState"] = ("==", self._entryActive)
		collection = self._collections[collectionIndex]

		return reduce(
			np.logical_and,
			(
				self._queryOperations[operator](
					collection[attrName],
					queryValue
					)
				for attrName, (operator, queryValue) in operations.viewitems()
			)
		)


	def objectsByGlobalIndex(self, globalIndexes):
		"""Return a _UniqueObjectSet proxy for the objects (molecules) with the
		given global indexes (that is, global over all named collections).
		"""
		return _UniqueObjectSet(self, globalIndexes)


	def objectByGlobalIndex(self, globalIndex):
		"""Return a _UniqueObject proxy for the object (molecule) with the
		given global index (that is, global over all named collections).
		"""
		return _UniqueObject(self, globalIndex)


	def objectNames(self):
		"""Return a tuple of the object names (AKA molecule type names or
		collection names) sorted by name.
		"""
		return self._names

	def counts(self, collectionNames=None):
		"""
		Get the counts of objects for each collection name.

		Parameters:
			collectionNames (iterable of str): The collection names. If None
				(default), then all collections are counted in sorted-name order.

		Returns:
			A vector (1D numpy.ndarray) of counts in the order of the requested
			collectionNames or else in the order of `self.objectNames()`.
		"""
		object_counts = np.array(
			[(x["_entryState"] == self._entryActive).sum() for x in self._collections])

		if collectionNames is None:
			return object_counts
		else:
			return object_counts[self._collectionNamesToIndexes(collectionNames)]

	def _collectionNamesToIndexes(self, collectionNames):
		"""
		Convert an iterable of collection names into their corresponding
		indices into the ordered list of collections.

		Parameters:
			collectionNames (iterable of str): The collection names.

		Returns:
			An array of indices (non-negative integers).
		"""
		return np.array([self._nameToIndexMapping[name] for name in collectionNames])

	def _copy_specs(self):
		"""Return a copy of the collection specifications without bookkeeping
		specs, as suitable for constructing a new UniqueObjectsContainer.
		"""
		specifications = deepcopy(self._specifications)
		specs_to_remove = self._defaultSpecification.keys()
		for moleculeName, moleculeSpecs in specifications.iteritems():
			for spec in specs_to_remove:
				moleculeSpecs.pop(spec)
		return specifications

	def emptyLike(self):
		"""Return a new container with the same specs, akin to np.zeros_like()."""
		specifications = self._copy_specs()
		new_copy = UniqueObjectsContainer(specifications)
		return new_copy

	def __eq__(self, other):
		# TODO(jerry): Ignore inactive entries and index values.
		# TODO(jerry): Don't access other's private fields.
		if not isinstance(other, UniqueObjectsContainer):
			return False
		if self._specifications != other._specifications:
			return False
		for (selfCollection, otherCollection) in izip(self._collections, other._collections):
			if not np.array_equal(selfCollection, otherCollection):
				return False
		return True

	def __ne__(self, other):
		# assertNotEquals() calls `!=`.
		return not (self == other)


	def loadSnapshot(self, other):
		"""Load data from a snapshot UniqueObjectsContainer which must have the
		same specifications, otherwise there's been a schema change.
		"""
		assert isinstance(other, UniqueObjectsContainer)

		if self._specifications != other._specifications:
			raise ValueError(
				'Schema change loading a UniqueObjectsContainer snapshot',
				other._specifications, self._specifications)

		assert self._names == other._names  # double check

		self._globalReference = other._globalReference.copy()
		for index, collection in enumerate(other._collections):
			self._collections[index] = collection.copy()


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			collectionNames = self._names
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			**dict(
				zip(self._names, self._collections)
				+ [("_globalReference", self._globalReference)]
				) # TODO: consider caching this dict
			)


	def tableLoad(self, tableReader, tableIndex):
		for fieldName, value in tableReader.readRow(tableIndex).viewitems():
			if fieldName == "_globalReference":
				self._globalReference = value

			else:
				self._collections[self._names.index(fieldName)] = value


def copy_if_ndarray(object):
	"""Copy an ndarray object or return any other type of object as is.
	Prevent making a view instead of a copy.  # <-- TODO(jerry): Explain.
	"""
	return object.copy() if isinstance(object, np.ndarray) else object


class _UniqueObject(object):
	"""
	A proxy for a row in a container, referring to a specific unique object.
	Primarily used as a way to manipulate individual molecules with a python
	object-like interface.
	"""

	__slots__ = ("_container", "_globalIndex", "_collectionIndex", "_objectIndex")


	def __init__(self, container, globalIndex):
		"""Construct a _UniqueObject proxy for the unique object (molecule) in
		the given container with the given global index.
		"""
		self._container = container
		self._globalIndex = globalIndex
		globalReference = container._globalReference[globalIndex]
		self._collectionIndex = globalReference["_collectionIndex"]
		self._objectIndex = globalReference["_objectIndex"]


	# def uniqueId(self):
	# 	return self.attr("_uniqueId")


	def attr(self, attribute):
		"""Return the named attribute of the unique object."""
		entry = self._container._collections[self._collectionIndex][self._objectIndex]

		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")

		return copy_if_ndarray(entry[attribute])


	def attrs(self, *attributes):
		"""Return a tuple containing the named attributes of the unique object."""
		entry = self._container._collections[self._collectionIndex][self._objectIndex]

		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")

		return tuple(copy_if_ndarray(entry[attribute]) for attribute in attributes)


	def attrIs(self, **attributes):
		"""Set named attributes of the unique object."""
		entry = self._container._collections[self._collectionIndex][self._objectIndex]

		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")

		for attribute, value in attributes.viewitems():
			if isinstance(entry[attribute], np.ndarray):
				# Fix for the circumstance that the attribute is an ndarray -
				# without the [:] assignment, only the first value will be
				# assigned (probably a NumPy bug)
				entry[attribute][:] = value

			else:
				entry[attribute] = value


	def __hash__(self):
		# TODO: replace with unique ID
		return hash((self._collectionIndex, self._objectIndex))


	def __eq__(self, other):
		if not isinstance(other, _UniqueObject):
			return False

		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return self._globalIndex == other._globalIndex


	def __ne__(self, other):
		return not self.__eq__(other)


class _UniqueObjectSet(object):
	"""
	An ordered sequence of _UniqueObject proxies for unique objects (molecules).
	Iterable and ordered.  Accessors allow for manipulating sets in lump.
	Internally this stores the objects' global indexes.
	"""

	def __init__(self, container, globalIndexes):
		"""Construct a _UniqueObjectSet for unique objects (molecules) in the
		given container with the given global indexes. The result is an
		iterable, ordered sequence (not really a set).
		"""
		self._container = container
		self._globalIndexes = np.array(globalIndexes, np.int)


	def __contains__(self, uniqueObject):
		if not self._container is uniqueObject._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return uniqueObject.attr("_globalIndex") in self._globalIndexes


	def __iter__(self):
		return (_UniqueObject(self._container, globalIndex)
			for globalIndex in self._globalIndexes)


	def __len__(self):
		return self._globalIndexes.size


	def __eq__(self, other):
		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return np.array_equal(self._globalIndexes, other._globalIndexes)

	def __ne__(self, other):
		return not self.__eq__(other)


	# TODO: More set-like operations (intersection, etc.)
	def __or__(self, other):
		# Use `|` -- not `or` -- to call this set operation.
		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return _UniqueObjectSet(
			self._container,
			np.lib.arraysetops.union1d(self._globalIndexes, other._globalIndexes)
			)


	def __getitem__(self, index):
		return _UniqueObject(self._container, self._globalIndexes[index])


	# def uniqueIds(self):
	# 	return self.attr("_uniqueId")


	def attr(self, attribute):
		"""Return the named attribute for all objects in this collection."""
		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		container = self._container
		globalReference = container._globalReference
		if (globalReference["_entryState"][self._globalIndexes] == container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		attributeDtype = container._collections[uniqueColIndexes[0]].dtype[attribute]

		values = np.zeros(
			self._globalIndexes.size,
			dtype = attributeDtype
			)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			values[globalObjIndexes] = container._collections[collectionIndex][attribute][objectIndexesInCollection]

		return values


	def attrs(self, *attributes):
		"""Return a tuple containing the named attributes for all objects in
		this sequence.
		"""
		return tuple(
			self.attr(attribute) for attribute in attributes
			)


	def attrsAsStructArray(self, *attributes):
		"""Return a structured array containing the named (or all) attributes
		from all objects in this sequence.
		"""
		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		container = self._container
		globalReference = container._globalReference
		if (globalReference["_entryState"][self._globalIndexes] == container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		if len(attributes) == 0:
			attributes = container._collections[uniqueColIndexes[0]].dtype.names

		attributes = list(attributes)

		attributeDtypes = [
			(attribute, container._collections[uniqueColIndexes[0]].dtype[attribute])
			for attribute in attributes
			]

		values = np.zeros(
			self._globalIndexes.size,
			dtype = attributeDtypes
			)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			values[globalObjIndexes] = container._collections[collectionIndex][attributes][objectIndexesInCollection]

		return values

	def attrIs(self, **attributes):
		"""Set named attributes of all the unique objects in this sequence."""
		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		container = self._container
		globalReference = container._globalReference
		if (globalReference["_entryState"][self._globalIndexes] == container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			for attribute, values in attributes.viewitems():
				valuesAsArray = np.array(values, ndmin = 1)

				if valuesAsArray.shape[0] == 1: # is a singleton
					container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray

				else:
					container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray[globalObjIndexes]


	def delByIndexes(self, indexes):
		"""Delete unique objects by indexes into this sequence."""
		# TODO(jerry): This could just call a delete method on all its
		# _UniqueObject instances, which in turn could call objectDel() on the
		# container or split the work so they each update their private state.

		globalIndexes = self._globalIndexes[indexes]

		# TODO: cache these properties? should be static
		container = self._container
		globalReference = container._globalReference
		collectionIndexes = globalReference["_collectionIndex"][globalIndexes]
		objectIndexes = globalReference["_objectIndex"][globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			container._collections[collectionIndex][objectIndexesInCollection] = np.zeros(
				1, dtype=container._collections[collectionIndex].dtype)

		globalReference[globalIndexes] = np.zeros(1, dtype=globalReference.dtype)


def _partition(objectRequestsArray, requestNumberVector, requestProcessArray, randomState):
	# Arguments:
	# objectRequestsArray: 2D bool array, (molecule)x(request)
	# requestNumberVector: number of molecules request, by request
	# requestProcessArray: 2D bool array, (request)x(process)
	# Returns:
	# partitionedMolecules: 2D bool array, (molecule)x(process)

	# TODO: full documentation/writeup, better docstring

	# Build matrix for optimization

	nObjects = objectRequestsArray.shape[0]
	nRequests = requestNumberVector.size
	nProcesses = requestProcessArray.shape[1]

	if nProcesses == 0:
		# Return nothing
		return np.zeros((nObjects, nProcesses), np.bool)

	# Make into structured array to condense the problem into unique rows
	objectRequestsStructured = objectRequestsArray.view(
		dtype = objectRequestsArray.dtype.descr * nRequests)

	uniqueEntriesStructured, mapping = np.unique(objectRequestsStructured,
		return_inverse = True)

	uniqueEntries = uniqueEntriesStructured.view((np.bool, nRequests))

	counts = np.bincount(mapping) # the number of each condensed molecule type

	nObjectTypes = counts.size

	# Some index mapping voodoo
	where0, where1 = np.where(uniqueEntries)

	nConnections = where0.size

	argsort = np.argsort(where1)

	moleculeToRequestConnections = np.zeros((nObjectTypes + nRequests,
		nConnections), np.int64)

	upperIndices = (where0, np.arange(where1.size)[argsort])
	lowerIndices = (nObjectTypes + where1[argsort], np.arange(where1.size))
	# End voodoo

	moleculeToRequestConnections[upperIndices] = -1
	moleculeToRequestConnections[lowerIndices] = 1

	# Create the matrix and fill in the values
	matrix = np.zeros(
		(nObjectTypes + nRequests + nProcesses,
			nObjectTypes + nConnections + 2*nProcesses),
		np.int64
		)

	# Molecule "boundary fluxes"
	matrix[:nObjectTypes, :nObjectTypes] = np.identity(nObjectTypes)

	# Flow from molecule type to request
	matrix[:nObjectTypes + nRequests,
		nObjectTypes:nObjectTypes+nConnections] = moleculeToRequestConnections

	# Flow from request to process
	matrix[nObjectTypes:nObjectTypes+nRequests,
		nObjectTypes+nConnections:nObjectTypes+nConnections+nProcesses][np.where(requestProcessArray)] = -requestNumberVector

	matrix[nObjectTypes + nRequests:,
		nObjectTypes+nConnections:nObjectTypes+nConnections+nProcesses] = np.identity(nProcesses)

	# Process "boundary fluxes"
	matrix[nObjectTypes + nRequests:,
		-nProcesses:] = -np.identity(nProcesses)

	# Create other linear programming parameters

	objective = np.zeros(matrix.shape[1], np.float)
	objective[-nProcesses:] = 1 # objective is to maximize process satisfaction
	# TODO: experiment with non-unity process weightings

	b = np.zeros(matrix.shape[0], np.float) # conservation law, i.e. b = 0 = Ax

	lowerBound = np.zeros(matrix.shape[1], np.float) # matrix is defined such that all values are >= 0

	upperBound = np.empty(matrix.shape[1], np.float)
	upperBound[:] = np.inf
	upperBound[:nObjectTypes] = counts # can use up to the total number of molecules
	upperBound[-nProcesses:] = 1 # processes can be up to 100% satisfied

	# Optimize

	solution = lp.linearProgramming(
		"maximize", objective,
		matrix.astype(np.float), b, # cvxopt requres floats
		lowerBound, upperBound,
		"S", "C", # no idea what these are supposed to do
		None # no options
		)[0].flatten()

	# Convert solution to amounts allocated to each process

	countsOfMoleculeTypeByConnection = -moleculeToRequestConnections[:nObjectTypes, :] * solution[nObjectTypes:nObjectTypes+nConnections]

	countsOfMoleculeTypeByConnection = np.floor(countsOfMoleculeTypeByConnection) # Round down to prevent oversampling

	countsOfMoleculeTypeByRequests = np.dot(
		countsOfMoleculeTypeByConnection,
		moleculeToRequestConnections[nObjectTypes:, :].T
		)

	countsOfMoleculeTypeByProcess = np.dot(
		countsOfMoleculeTypeByRequests,
		requestProcessArray
		)

	processOffsetsOfSelectedObjects = np.c_[
		np.zeros(nObjectTypes),
		np.cumsum(countsOfMoleculeTypeByProcess, axis = 1)
		].astype(np.int64)

	# TODO: find a way to eliminate the for-loops!
	partitionedMolecules = np.zeros((nObjects, nProcesses), np.bool)

	for moleculeIndex in np.arange(uniqueEntriesStructured.size):
		indexesOfSelectedObjects = np.where(moleculeIndex == mapping)[0]
		randomState.shuffle(indexesOfSelectedObjects)

		for processIndex in np.arange(nProcesses):
			start = processOffsetsOfSelectedObjects[moleculeIndex, processIndex]
			stop = processOffsetsOfSelectedObjects[moleculeIndex, processIndex + 1]
			selectedIndexes = indexesOfSelectedObjects[start : stop]

			partitionedMolecules[
				selectedIndexes,
				processIndex] = True

	return partitionedMolecules
