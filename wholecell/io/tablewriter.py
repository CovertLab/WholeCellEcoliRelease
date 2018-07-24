
from __future__ import absolute_import
from __future__ import division

import os
import json

import numpy as np

from wholecell.utils import filepath

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
		filepath.makedirs(path)

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
	"""
	Generic live output writer for NumPy ndarrays.

	NumPy can save and load arrays to disk.  This class provides a convenient
	interface to repeated appending of NumPy array data to an output file,
	which can be loaded one entry at a time or all at once via the companion
	class, TableReader.

	Output file structure:

	<root directory> : Root path, provided by during instantiation.
		/attributes : Directory for saved JSON data, provided by the user.
			/<attribute name> : A JSON file.
		/columns : Directory for the main output streams ("columns").
			/<column name> : Directory for a specific column.
				data : Contains the JSON-serialized format specification for
					the Numpy array dtype and appended binary data from NumPy
					ndarrays.
				offsets : A list of integers (one per line) corresponding to
					the associated byte offset for each entry in the column.
		/metadata : Directory for saving any critical internals for TableWriter
				and TableReader.
			version : An integer representing the "version number" of the data
				format specification.

	Parameters
	----------
	path : str
		Path to the output location (a directory) where data will be saved.  If
			 he directory already exists, no error will be thrown.

	See also
	--------
	whoelcell.io.tablereader.TableReader

	Notes
	-----

	Data written to columns can be of fixed or variable size.  If fixed, the
	output can be read in as a single, higher dimensional array by TableReader.
	Otherwise output will need to be handled one line at a time.

	0D and 1D array writing is supported.  Higher dimensions are outside of
	spec, and will likely be flattened when read by TableReader.
	TODO (John): throw an error for > 1D?

	Both simple and structured ndarray dtypes are supported.  Structured arrays
	are a good way to work around the dimensional limitations.

	NumPy object arrays (e.g. arrays not of one or more standard dtypes) are
	not supported.  They might work, but under the hood will rely on pickling
	for saving and loading, which will be terribly slow.

	TODO (John): Consider moving the Numpy dtype format specification out of
		the 'data' file and into its own file.  This was done originally to
		keep the format tightly associated with the data, improving portability
		and consistency with the np.save implementation.  However it leads to
		logical complications in saving and loading, and it mixes JSON string
		data with array binary data in one file.

	TODO (John): Drop the 'version' idea as well as the 'metadata' directory.
		We so rarely try to run new scripts on old data (or vice-versa) that
		I'm not sure this needs to exist or be maintained.

	TODO (John): Move the _Column class into the TableWriter namespace.
		There's no reason for it to be available at the module level.

	TODO (John): Instead of saving many files and sub-directories under one
		directory, save everything to an uncompressed archive.  Built-in module
		zipfile seems like a good option; np.savez uses it.

	TODO (John): Unit tests.

	"""

	def __init__(self, path):

		dirMetadata = filepath.makedirs(path, DIR_METADATA)

		open(os.path.join(dirMetadata, FILE_VERSION), "w").write(VERSION)

		self._dirAttributes = filepath.makedirs(path, DIR_ATTRIBUTES)
		self._attributeNames = []

		self._dirColumns = filepath.makedirs(path, DIR_COLUMNS)
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
