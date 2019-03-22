
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
	"""
	Base exception class for TableWriter-associated exceptions.
	"""
	pass

class MissingFieldError(TableWriterError):
	"""
	An error raised when TableWriter.append is called without providing all
	field names.
	"""
	pass

class UnrecognizedFieldError(TableWriterError):
	"""
	An error raised when TableWriter.append is called with an unexpected field
	name.
	"""
	pass

class AttributeAlreadyExistsError(TableWriterError):
	"""
	An error raised when TableWriter.writeAttributes is called with an
	attribute name that has already been used.
	"""
	pass

class AttributeTypeError(TableWriterError):
	"""
	An error raised when TableWriter.writeAttributes is called with a type that
	does not appear to be JSON-serializable.
	"""
	pass


class _Column(object):
	"""
	Manages the written data for a specific TableWriter field.

	Each field in a 'table' corresponds to a 'column' that is written to on
	each append operation.  This private class encapsulates the logic and data
	for a particular 'column'.

	Parameters
	----------
	path : str
		The path to the sub-directory associated with this particular column.

	Notes
	-----
	See TableWriter for more information about output file and directory
	structure.

	TODO (John): With some adjustment this class could be made public, and used
		as a lightweight alternative to TableWriter in addition to part of
		TableWriter's internal implementation.

	"""

	def __init__(self, path):
		filepath.makedirs(path)

		self._data = open(os.path.join(path, FILE_DATA), "w")
		self._offsets = open(os.path.join(path, FILE_OFFSETS), "w")

		self._dtype = None


	def append(self, value):
		"""
		Appends an array-like to the end of a column.

		On first call, the NumPy dtype of the column is inferred from the
		input.  On subsequent calls this dtype is checked for consistency.

		Parameters
		----------
		value : array-like
			A NumPy ndarray or anything that can be cast to an array via
			np.asarray, including scalars (i.e. 0D arrays).

		"""

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
		"""
		Close the files associated with the column.

		While running, each column keeps two files open; one associated with
		the data itself (as well as the dtype information), and one associated
		with the offsets between entries.  This method explicitly closes those
		files.

		Notes
		-----
		Trying to append after closing will raise an error.

		"""

		self._data.close()
		self._offsets.close()


	def __del__(self):
		"""
		Explicitly closes the output files once the instance is totally
		dereferenced.

		Notes
		-----
		This will lead to errors consequent of operating on a closed file if
		references to the output files (which ought to be private attributes)
		exist.  This is desirable, because such references should not be made.

		TODO (John): Decide whether we want this sort of 'irresponsible
			developer' error prevention.  E.g. see the Wikipedia page on
			"offensive programming".

		"""
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
		Path to the directory that will be created.  All data will be saved
		under this directory.  If the directory already exists, no error will
		occur.

	See also
	--------
	wholecell.io.tablereader.TableReader

	Notes
	-----
	The terms used in this class are adapted from the analogy of a spreadsheet
	or table.  Each append operation adds a new 'row' to the table, also called
	an 'entry'.  Every 'column' corresponds to a 'field' and contains a
	particular set of data (i.e. of a fixed type) for all entries.

	Attributes are meant for user-provided annotation that may be useful for
	downstream analysis or otherwise improve portability.  E.g. a list of names
	associated with the elements of a vector-column.

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

	TODO (John): Test portability across machines (particularly, different
	operating systems).

	TODO (John): Consider separating out the fixed and variable size
		implementations.  Further, consider writing all fields simultaneously
		as part of a structured array (i.e. a hybrid data type).

	"""

	def __init__(self, path):

		dirMetadata = filepath.makedirs(path, DIR_METADATA)

		open(os.path.join(dirMetadata, FILE_VERSION), "w").write(VERSION)

		self._dirAttributes = filepath.makedirs(path, DIR_ATTRIBUTES)
		self._attributeNames = []

		self._dirColumns = filepath.makedirs(path, DIR_COLUMNS)
		self._columns = None


	def append(self, **namesAndValues):
		"""
		Write a new set of values to each column.

		On the first call to this method, the columns will be set up with the
		appropriate names and data type information.  Subsequent calls will
		validate that the names and data types are consistent.

		Parameters
		----------
		**namesAndValues : dict of {string: array-like} pairs
			The column names (fields) and associated values to append to the
			end of the columns.

		Notes
		-----
		All fields must be provided every time this method is called.

		"""

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
		"""
		Writes JSON-serializable data.

		Oftentimes some additional data is needed to contextualize the data in
		the columns.  This method can be called to write JSON-serializable
		data (e.g. a list of strings) alongside the column data.

		Parameters
		----------
		**namesAndValues : dict of {string: JSON-serializable} pairs
			The attribute names and associated values.

		Notes
		-----
		This method can be called at any time, so long as an attribute name is
		not reused.

		"""

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
		"""
		Close the output files (columns).

		Notes
		-----
		Trying to append after closing will raise an error.

		"""

		if self._columns is not None:
			for column in self._columns.viewvalues():
				column.close()


	def __del__(self):
		"""
		Explicitly closes the output files once the instance is totally
		dereferenced.

		Notes
		-----
		TODO (John): Decide whether we want to implement __del__ (see analogous
			method in the _Column class for more details).

		"""
		self.close()
