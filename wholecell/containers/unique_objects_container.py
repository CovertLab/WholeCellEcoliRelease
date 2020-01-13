"""
A UniqueObjectsContainer tracks the attributes of unique objects, which are
typically molecules.

The UniqueObjectsContainer uses _UniqueObject instances to present O-O proxies
to specific molecules.

SEE THE decomp2() CAUTION ABOUT PERSISTENT DATA (PICKLE) FORMAT.
"""

from __future__ import absolute_import, division, print_function

from copy import deepcopy
from itertools import izip, product
from enum import Enum

import numpy as np
import zlib

# TODO: object transfer between UniqueObjectsContainer instances

ZLIB_LEVEL = 7

Access = Enum('Access', 'EDIT DELETE')

class UniqueObjectsContainerException(Exception):
	pass

class UniqueObjectsPermissionException(Exception):
	pass

class UniqueObjectsMergeConflictException(Exception):
	pass

class UniqueObjectsInvalidSubmassNameException(Exception):
	pass


def decomp(specifications, compressed_collections, global_ref_count, submass_diff_names=None):
	raise UniqueObjectsContainerException("Can't read an old format container")

def decomp2(specifications, compressed_collections, global_ref_count,
		next_unique_index, submass_diff_names=None):
	"""Decompress the arguments into a UniqueObjectsContainer. The name is
	intentionally short and awkward for intended use limited to pickling. It
	calls the constructor to set up indexes and caches, unlike `__setstate__`.

	CAUTION: Changes to this class must manage forward and backward
	compatibility of persistent data via decomp and __reduce__.
	* For new code reading old data, make decomp either properly read and
	upgrade old formats or fail fast on old formats.
	* For old code reading new data, make __reduce__ format the new data so
	old code will either read it properly or fail fast.
	* Changes to instance structure also affect __eq__ and loadSnapshot.

	The current stored format is sensitive to all the internal representation
	details including inactive entries, bookkeeping fields, and index fields
	which could all be reconstructed here with a bunch of work and a more
	complex __eq__() method.

	Args:
		specifications (dict): dtype specs as passed to UniqueObjectsContainer()
		compressed_collections (list[bytes]): zlib-compressed bytes of the collections ndarrays
		global_ref_count (int): the size of the global references ndarray
		next_unique_index (int): the container's next globally unique index
		submass_diff_names (optional, list[strings]): list of the names of
			submass difference attributes

	Returns:
		A filled-in UniqueObjectsContainer.
	"""
	container = UniqueObjectsContainer(
		specifications, submass_diff_names=submass_diff_names)

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

	container._next_unique_index = next_unique_index
	return container


def compress_ndarray(array):
	return zlib.compress(array.tobytes(), ZLIB_LEVEL)

def decompress_ndarray(compressed_bytes, dtype):
	return np.frombuffer(zlib.decompress(compressed_bytes), dtype).copy()


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

	# Bookkeeping fields to add to every struct type
	_defaultSpecification = {
		"_entryState": np.int8,  # See state descriptions above
		"_globalIndex": np.int64,  # Index of object in the _globalReference array (collection)
		"unique_index": np.int64,  # Unique index assigned to each object
		}

	_globalReferenceDtype = {
		"_entryState":np.int8,  # See state descriptions above
		"_collectionIndex":np.int64,  # Index of collection in self._collections
		"_objectIndex":np.int64,  # Index of object in each structured array
		}

	_fractionExtendEntries = 0.1 # fractional rate to increase number of entries in the structured array (collection)


	def __init__(self, specifications, submass_diff_names=None):
		"""
		Parameters:
			specifications (Dict[str, Dict[str, str]]): Maps the unique molecule
				names (collection names) to the {attribute_name: dtype} molecule
				attributes (structured array fields). The dtype declarations are
				strings describing NumPy scalar types.

				Example: {
					'DNA polymerase': {'bound': 'bool', 'location': 'int32'},
					'RNA polymerase': {'bound': 'bool', 'location': 'int32', 'decay': 'float32'},
					}

			submass_diff_names (optional, List[str]): List of attribute names
				that correspond to added masses of the unique molecule.
				TODO (ggsun): Ideally, the container should be agnostic to
					this information. This list is passed to distinguish
					edits to the container that would change the mass of the
					cell.
		"""
		self._collections = [] # ordered list of numpy structured arrays
		self._nameToIndexMapping = {} # collectionName:index of associated structured array

		self._specifications = deepcopy(specifications) # collectionName:{attributeName:type}
		self._names = tuple(sorted(self._specifications.keys())) # sorted collection names

		self.submass_diff_names_list = submass_diff_names or []
		self.submass_diff_names_set = frozenset(self.submass_diff_names_list)

		# List of requests
		self._requests = []

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
		# _collections, i.e. an index across all the unique object names
		# (molecule types) in this container.
		self._globalReference = np.zeros(
			1, dtype = make_dtype_spec(self._globalReferenceDtype))

		# Initialize next unique index to zero
		self._next_unique_index = 0


	def __reduce__(self):
		"""
		Reduce the container to its defining state for pickling.
		Compress the state for transmission efficiency.
		Return a callable object and its args.
		"""
		# TODO(jerry): Compress the specs?
		# TODO(jerry): Squeeze out inactive entries and _entryState fields, but
		#   indicate which slots are inactive or else make __eq__ ignore
		#   differences in their positions and impact on indexes. Worth the work?
		if len(self._requests) != 0:
			raise UniqueObjectsContainerException(
				"Cannot pickle container with unapplied requests. Run .merge() to apply the requests before pickling."
				)

		specs = self._copy_specs()
		compressed_collections = [compress_ndarray(col) for col in self._collections]
		return decomp2, (
			specs, compressed_collections, self._globalReference.size,
			self._next_unique_index, self.submass_diff_names_list)


	def __eq__(self, other):
		# TODO(jerry): Ignore inactive entries and index values.
		if not isinstance(other, UniqueObjectsContainer):
			return False
		if self._next_unique_index != other._next_unique_index:
			return False
		if self._specifications != other._specifications:
			return False
		if self.submass_diff_names_list != other.submass_diff_names_list:
			return False
		for (selfCollection, otherCollection) in izip(self._collections, other._collections):
			if not np.array_equal(selfCollection, otherCollection):
				return False
		return True


	def __ne__(self, other):
		# assertNotEquals() calls `!=`.
		return not (self == other)


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
		"""
		Add nObjects new objects/molecules of the named type, all with the
		given attributes. The objects added here do not wait for merge and are
		added immediately. Return the unique indexes of the newly added
		molecules.
		"""
		# Add unique index to dict of attributes
		unique_indexes = np.arange(
			self._next_unique_index, self._next_unique_index + nObjects)
		attributes['unique_index'] = unique_indexes

		# Increment value for next available unique index
		self._next_unique_index += nObjects

		# Add molecules and return unique indexes
		self._add_new_objects(collectionName, nObjects, attributes)
		return unique_indexes


	def objectNew(self, collectionName, **attributes):
		"""
		Add a new object/molecule of the named type with the given attributes.
		"""
		self.objectsNew(collectionName, 1, **attributes)


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


	def objects(self, access=()):
		"""
		Return a _UniqueObjectSet proxy for all active objects (molecules).
		The access argument determines the level of access permission the proxy
		has to the molecules in the container.
		"""
		return _UniqueObjectSet(self,
			np.where(self._globalReference["_entryState"] == self._entryActive)[0],
			access=access
			)


	def objectsInCollection(self, collectionName, process_index=None,
			access=()):
		"""
		Return a _UniqueObjectSet proxy for all objects (molecules) belonging
		to a named collection. The access argument determines the level of
		access permission the proxy has to the molecules in the container.
		"""
		collectionIndex = self._nameToIndexMapping[collectionName]

		active_mask = self._find_active_entries(collectionIndex)

		return _UniqueObjectSet(self,
			self._collections[collectionIndex]["_globalIndex"][active_mask],
			process_index=process_index,
			access=access
			)


	def objectsInCollections(self, collectionNames, process_index=None,
			access=()):
		"""
		Return a _UniqueObjectSet proxy for all objects (molecules) belonging
		to one of the given collection names. The access argument determines
		the level of access permission the proxy has to the molecules in the
		container.
		"""
		collectionIndexes = [self._nameToIndexMapping[collectionName] for collectionName in collectionNames]
		active_masks = []

		for collectionIndex in collectionIndexes:
			active_masks.append(self._find_active_entries(collectionIndex))

		return _UniqueObjectSet(self,
			np.concatenate([
			self._collections[collectionIndex]["_globalIndex"][result]
			for collectionIndex, result in izip(collectionIndexes, active_masks)]),
			process_index=process_index,
			access=access
			)


	def _find_active_entries(self, collectionIndex):
		"""
		For the structured array at the given index, return a boolean array
		that's True where entries are active.
		"""
		collection = self._collections[collectionIndex]
		return collection["_entryState"] == self._entryActive


	def objectsByGlobalIndex(self, globalIndexes, access=()):
		"""Return a _UniqueObjectSet proxy for the objects (molecules) with the
		given global indexes (that is, global over all named collections).
		"""
		return _UniqueObjectSet(self, globalIndexes, access=access)


	def objectByGlobalIndex(self, globalIndex, access=()):
		"""Return a _UniqueObject proxy for the object (molecule) with the
		given global index (that is, global over all named collections).
		"""
		return _UniqueObject(self, globalIndex, access=access)


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
		new_copy = UniqueObjectsContainer(specifications, self.submass_diff_names_list)
		return new_copy


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

		self._next_unique_index = other._next_unique_index


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


	def add_request(self, **fields):
		"""
		Adds a request made from a _UniqueObjectSet instance to the list of
		requests to handle. fields["type"] can be "edit", "submass", "delete",
		or "new_molecule". If the request type is "new_molecule", adds unique
		indexes to the attributes of the new molecules, and returns the
		indexes.
		"""
		if fields['type'] == 'new_molecule':
			# Add unique index to dict of attributes
			unique_indexes = np.arange(
				self._next_unique_index,
				self._next_unique_index + fields['nObjects'])
			fields['attributes']['unique_index'] = unique_indexes

			# Increment value for next available unique index
			self._next_unique_index += fields['nObjects']

			# Append request and return unique indexes
			self._requests.append(fields)
			return unique_indexes

		else:
			# Append request and return None
			self._requests.append(fields)
			return None


	def _add_new_objects(self, collectionName, nObjects, attributes):
		"""
		Adds new objects to array with the given initial attributes. All new
		objects are given unique indexes.
		"""
		collectionIndex = self._nameToIndexMapping[collectionName]
		objectIndexes, globalIndexes = self._getFreeIndexes(collectionIndex, nObjects)

		collection = self._collections[collectionIndex]

		# TODO(jerry): Would it be faster to copy one new entry to all rows
		# then set the Index columns?
		collection["_entryState"][objectIndexes] = self._entryActive
		collection["_globalIndex"][objectIndexes] = globalIndexes

		for attrName, attrValue in attributes.viewitems():
			collection[attrName][objectIndexes] = attrValue

		self._globalReference["_entryState"][globalIndexes] = self._entryActive
		self._globalReference["_collectionIndex"][globalIndexes] = collectionIndex
		self._globalReference["_objectIndex"][globalIndexes] = objectIndexes


	def _delete_objects(self, global_indexes):
		"""
		Deletes objects referenced by global_indexes. Returns the collection
		indexes of the deleted molecules, and the sum of submass differences
		of all the delete molecules.
		"""
		self._check_deleted_objects(global_indexes)

		global_reference = self._globalReference

		collection_indexes = global_reference[global_indexes]["_collectionIndex"]
		object_indexes = global_reference[global_indexes]["_objectIndex"]

		deleted_submasses = np.zeros(len(self.submass_diff_names_list))

		unique_col_indexes, inverse = np.unique(collection_indexes,
			return_inverse=True)

		for idx, collection_index in enumerate(unique_col_indexes):
			globalObjIndexes = np.where(inverse == idx)
			objectIndexesInCollection = object_indexes[globalObjIndexes]

			for submass_name_idx, submass_diff_name in enumerate(
					self.submass_diff_names_list):
				deleted_submasses[submass_name_idx] += self._collections[
					collection_index][submass_diff_name][
					objectIndexesInCollection].sum()

			self._collections[collection_index][
				objectIndexesInCollection] = np.zeros(
				1, dtype=self._collections[collection_index].dtype)

		global_reference[global_indexes] = np.zeros(1,
			dtype=global_reference.dtype)

		return collection_indexes, deleted_submasses


	def update_attribute(self, global_indexes, attributes, scale):
		"""
		Updates the values of the existing attributes of objects. If scale=0,
		the values are overwritten by the values given in the attributes
		argument. If scale=1, the values given as arguments are added to the
		existing values.
		"""
		self._check_deleted_objects(global_indexes)

		globalReference = self._globalReference

		collectionIndexes = globalReference["_collectionIndex"][global_indexes]
		objectIndexes = globalReference["_objectIndex"][global_indexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes,
			return_inverse=True)

		for idx, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == idx)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			for attribute, deltas in attributes.viewitems():
				deltas_as_array = np.array(deltas, ndmin=1)
				values = self._collections[collectionIndex][attribute][
					objectIndexesInCollection]

				if deltas_as_array.shape[0] == 1:  # is a singleton
					self._collections[collectionIndex][attribute][
						objectIndexesInCollection] = values*scale + deltas_as_array

				else:
					self._collections[collectionIndex][attribute][
						objectIndexesInCollection] = values*scale + deltas_as_array[globalObjIndexes]


	def _check_deleted_objects(self, global_indexes):
		"""
		Checks if any of the objects referenced by the given global indexes
		have already been removed. Raises exception if this is the case.
		"""
		globalReference = self._globalReference

		if (globalReference["_entryState"][
				global_indexes] == self._entryInactive).any():
			raise UniqueObjectsContainerException(
				"One or more object was deleted from the set")



	def merge(self):
		"""
		Loops through the list of all requests and makes the requested changes.
		Raises exception if conflicting requests are made on the same object.
		Returns a copied list of requests to the UniqueMolecules state, which
		then is used to compute the mass changes that occurred in each process.
		The list of requests is emptied after the copy is made.
		"""
		resolver = []

		priority = ["new_molecule", "edit", "submass", "delete"]
		sorted_requests = sorted(self._requests, key=lambda k: priority.index(k["type"]))

		# Loop through all requests
		for req in sorted_requests:
			# Apply requested attribute edits
			if req["type"] == "edit":
				self.update_attribute(req["globalIndexes"], req["attributes"], 0)
				resolver.extend(
					list(product(req["globalIndexes"], req["attributes"].keys()))
					)

			# Apply requested submass edits
			elif req["type"] == "submass":
				self.update_attribute(req["globalIndexes"], req["added_masses"], 1)


			# Apply requested deletions
			elif req["type"] == "delete":
				collection_indexes, deleted_submasses = self._delete_objects(
					req["globalIndexes"]
					)

				# Add new keys to request (used to calculate mass difference)
				req["collection_indexes"] = collection_indexes
				req["deleted_submasses"] = deleted_submasses

			# Apply requested new molecule additions
			elif req["type"] == "new_molecule":
				self._add_new_objects(
					req["collectionName"], req["nObjects"], req["attributes"]
					)

		# Check for multiple edit requests on the same attribute of a same molecule
		if len(resolver) != len(set(resolver)):
			raise UniqueObjectsMergeConflictException(
				"Merge conflict detected - two processes attempted to edit same attribute of same unique molecule."
				)

		# Empty the original list
		self._requests = []

		return sorted_requests


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

	__slots__ = ("_container", "_globalIndex", "_collectionIndex", "_objectIndex", "_access")


	def __init__(self, container, globalIndex, access=()):
		"""Construct a _UniqueObject proxy for the unique object (molecule) in
		the given container with the given global index.
		"""
		self._container = container
		self._globalIndex = globalIndex
		globalReference = container._globalReference[globalIndex]
		self._collectionIndex = globalReference["_collectionIndex"]
		self._objectIndex = globalReference["_objectIndex"]
		self._access = access


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
		if Access.EDIT not in self._access:
			raise UniqueObjectsPermissionException(
				"Can't edit attributes of read-only objects."
			)

		# Submass attributes must be edited through specialized methods.
		if not self._container.submass_diff_names_set.isdisjoint(attributes.keys()):
			raise UniqueObjectsPermissionException(
				"Can't modify submass differences with attrIs(). Use add_submass_by_name() or add_submass_by_array() instead."
				)

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

	def __init__(self, container, globalIndexes, process_index=None,
			access=()):
		"""
		Construct a _UniqueObjectSet for unique objects (molecules) in the
		given container with the given global indexes. The result is an
		iterable, ordered sequence (not really a set). The access argument
		determines the level of access permission this instance has to the
		molecules in the container, with an empty tuple meaning that the
		instance has read-only access. If an instance is initialized from a
		process View, self._process_index is set to the index of the process.
		In this case, attrIs(), add_submass_by_name(), and
		add_submass_by_array() all submit an edit request to the container, and
		the container waits until merge to actually perform the edits. If
		self._process_index is not set, the edits are made directly on the
		container.
		"""
		self._container = container
		self._globalIndexes = np.array(globalIndexes, np.int)
		self._process_index = process_index
		self._access = access


	def __contains__(self, uniqueObject):
		if not self._container is uniqueObject._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return uniqueObject.attr("_globalIndex") in self._globalIndexes


	def __iter__(self):
		return (_UniqueObject(self._container, globalIndex, access=self._access)
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
		return _UniqueObject(self._container, self._globalIndexes[index],
			access=self._access)


	def set_access_level(self, access=()):
		"""
		Resets the access permission to the underlying container given to the
		UniqueObjectSet instance.
		"""
		self._access = access


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

		for idx, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == idx)
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

		for idx, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == idx)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			values[globalObjIndexes] = container._collections[collectionIndex][attributes][objectIndexesInCollection]

		return values


	def attrIs(self, **attributes):
		"""
		Set named attributes of all the unique objects in this sequence.
		This is not permitted for read-only sets.
		"""
		if Access.EDIT not in self._access:
			raise UniqueObjectsPermissionException(
				"Can't edit attributes of read-only objects."
			)

		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		# Submass attributes must be edited through specialized methods.
		if not self._container.submass_diff_names_set.isdisjoint(attributes.keys()):
			raise UniqueObjectsPermissionException(
				"Can't modify submass differences with attrIs(). Use add_submass_by_name() or add_submass_by_array() instead."
				)

		# Submit edit request to container
		if self._process_index is not None:
			self._container.add_request(
				type="edit",
				globalIndexes=self._globalIndexes,
				process_index=self._process_index,
				attributes=attributes,
				)

		# Make edit directly on container now
		else:
			self._container.update_attribute(
				self._globalIndexes, attributes, 0
				)


	def add_submass_by_name(self, submass_name, delta_mass):
		"""
		Adds a given amount of mass to a specific submass type whose name is
		given as the argument. This should be used when a single specific
		submass needs to be added to unique molecules.
		- submass_name (str): name of submass being added (e.g. "DNA")
		- delta_mass (1D array, length equal to number of objects in set): mass
		being added to each unique object
		"""
		if Access.EDIT not in self._access:
			raise UniqueObjectsPermissionException(
				"Can't edit attributes of read-only objects."
				)

		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		submass_attr_name = "massDiff_" + submass_name

		if submass_attr_name not in self._container.submass_diff_names_set:
			raise UniqueObjectsInvalidSubmassNameException(
				'"%s" is not a valid submass name.' % (submass_name, )
			)

		added_masses = {
			submass_attr_name: delta_mass
			}

		# If the call comes from a process, submit request to container
		if self._process_index is not None:
			self._container.add_request(
				type="submass",
				globalIndexes=self._globalIndexes,
				process_index=self._process_index,
				added_masses=added_masses,
				)

		# Otherwise make edit directly on container
		else:
			self._container.update_attribute(
				self._globalIndexes, added_masses, 1
				)


	def add_submass_by_array(self, delta_mass):
		"""
		Adds a given amount of mass, given as an array, to all submass types.
		This should be used when multiple types of submasses need to be added
		to unique molecules simultaneously.
		- delta_mass (ndarray, with shape (N, M), N: number of objects in set
			M: number of submass types)
			: array of submasses being added to each object.
		"""
		if Access.EDIT not in self._access:
			raise UniqueObjectsPermissionException(
				"Can't edit attributes of read-only objects."
				)

		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		added_masses = {}
		for i, submass_attr_name in enumerate(self._container.submass_diff_names_list):
			added_masses[submass_attr_name] = delta_mass[:, i]

		# If the call comes from a process, submit request to container
		if self._process_index is not None:
			self._container.add_request(
				type="submass",
				globalIndexes=self._globalIndexes,
				process_index=self._process_index,
				added_masses=added_masses,
				)

		# Otherwise make edit directly on container
		else:
			self._container.update_attribute(
				self._globalIndexes, added_masses, 1
				)


	def delByIndexes(self, indexes):
		"""
		Submits a request to delete unique objects by indexes into this
		sequence. This is not permitted for read-only sets.
		"""
		if Access.DELETE not in self._access:
			raise UniqueObjectsPermissionException(
				"Can't delete molecules from this object without delete access."
			)

		# TODO(jerry): This could just call a delete method on all its
		# _UniqueObject instances, which in turn could call objectDel() on the
		# container or split the work so they each update their private state.
		globalIndexes = self._globalIndexes[indexes]

		self._container.add_request(
			type="delete",
			globalIndexes=globalIndexes,
			process_index=self._process_index,
			)
