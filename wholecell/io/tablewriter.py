
from __future__ import absolute_import, division, print_function

import os
import json
import numpy as np
import struct
from typing import Any, Dict, List, Optional, Set
import zlib

from wholecell.utils import filepath
import six


__all__ = [
	"TableWriter",
	"TableWriterError",
	"MissingFieldError",
	"UnrecognizedFieldError",
	"AttributeAlreadyExistsError",
	"AttributeTypeError",
	"VariableEntrySizeError",
	"TableExistsError",
	]


# --- Persistent values and structures ---
VERSION = 3  # should update this any time there is a spec-breaking change

FILE_ATTRIBUTES = "attributes.json"

# Chunk type and size.
CHUNK_HEADER = struct.Struct('>4s I')
COLUMN_CHUNK_TYPE = b'COLM'  # column file's header chunk
VARIABLE_COLUMN_CHUNK_TYPE = b'VCOL'  # variable-length column's header chunk
BLOCK_CHUNK_TYPE = b'BLOC'   # data block chunk
ROW_SIZE_CHUNK_TYPE = b'RWSZ'  # row size block chunk

# Datatype of row size chunks (contains number of array elements for each row)
ROW_SIZE_CHUNK_DTYPE = np.uint32

# Datatype of variable-length columns (Set to float to support the use of
# np.nan as filler values when reading)
VARIABLE_COLUMN_DATA_DTYPE = np.float64

# Column header struct. See the pack() calls for field details.
COLUMN_STRUCT = struct.Struct('>2I 2H')

# Variable-length column header struct.
VARIABLE_COLUMN_STRUCT = struct.Struct('H')

COMPRESSION_TYPE_NONE = 0
COMPRESSION_TYPE_ZLIB = 1

V2_DIR_COLUMNS = "columns"  # format v2's directory of column files
# ----------------------------------------


# zlib's default compression level is 6.
#
# Measuring zlib compression of Mass/time, BulkMolecules/atpRequested, and
# BulkMolecules/counts, higher levels yielded a fraction of a percent space
# savings at substantial time costs. For BulkMolecules/counts, level 7 adds
# ~40% time, level 8 takes 5x, and level 9 takes 15x time.
#
# See wholecell/tests/io/measure_zlib.py
ZLIB_LEVEL = 6

# Pack enough entries into each data block to total about BLOCK_BYTES_GOAL
# bytes before compression.
#
# Compressing multiple entries as a block saves considerable space for small
# and medium size entries. E.g. compressing an 8-byte Mass/time entry doubles
# the size, and adding a chunk header gets it to 3x, while compressing 512 of
# them into a block gets the entire chunk down to 68% of the input size.
#
# BulkMolecules/atpRequested has 96-byte entries which compress individually to
# 30% including header. It gets down to 10% after packing 16 entries together,
# 9% when packing 32 entries together, and about 7% when packing 512 entries
# together.
#
# Packing also saves compression time and presumably I/O time.
# TODO: Measure I/O time at different block sizes.
#
# The column file format isn't designed for random access but a reader could
# skip chunk to chunk without decompressing them to get to a desired block, and
# adding a chunk index would enable block-level random access.
#
# zlib supports incrementally compressing a bytestring at a time to save buffers
# and use the cumulative data history to improve the overall compression ratio,
# but measurements show it can't obviate data blocks. For readers to decompress
# an entry at a time, the writer has to Z_SYNC_FLUSH each one which reduces the
# compression ratio (even if the writer discards the '\x00\x00\xff\xff' suffix
# and the reader restores it). So the writer still needs to group entries into
# blocks. Then cumulative incremental compression has a tiny +/- impact on
# compression size, and the ability to read such an entry without decompressing
# everything before it depends on saving & reloading the compression state,
# which takes space -- and that feature is not in the Python zlib library.
#
# See wholecell/tests/io/measure_zlib.py
BLOCK_BYTES_GOAL = 16384


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

class VariableEntrySizeError(TableWriterError):
	"""
	An error raised on attempt to write an entry that's not the same size as
	the previous entries in the same Column.
	"""
	pass

class TableExistsError(TableWriterError):
	"""
	An error raised on attempt to create a Table in a directory that already
	has a table (completed or in progress).
	"""
	pass


class _Column(object):
	"""
	Base class for managing the written data for a specific TableWriter field.

	Each field in a 'table' corresponds to a 'column' that is written to on
	each append operation.  This private class encapsulates the logic and data
	for a particular 'column'.

	The file format is a sequence of chunks. It begins with a COLUMN_CHUNK_TYPE
	or VARIABLE_COLUMN_CHUNK_TYPE header chunk followed by one or more
	BLOCK_CHUNK_TYPE data chunks (a column cannot be empty) containing
	optionally compressed NumPy ndarray data. In the case of variable-length
	columns, each data chunk is preceded by a row size chunk that specifies the
	sizes (number of array elements) for each row included in the data chunk.
	The format is extensible since readers should skip unrecognized chunk
	types.

	Parameters:
		path (str): The path for this particular column's data file.

	Notes
	-----
	See TableWriter for more information about output file and directory
	structure.  Each column creates a file containing a COLUMN_CHUNK_TYPE or
	a VARIABLE_COLUMN_CHUNK_TYPE chunk followed by the data BLOCK_CHUNK_TYPE
	chunks and ROW_SIZE_CHUNK_TYPE chunks. Each data block contains one or more
	NumPy array 'entries', optionally compressed. One can read chunks using the
	Python chunk library Chunk(file, False). A conforming reader must skip
	unrecognized chunks.

	TODO (John): With some adjustment this class could be made public, and used
		as a lightweight alternative to TableWriter in addition to part of
		TableWriter's internal implementation.
	"""
	def __init__(self, path, compression_type=COMPRESSION_TYPE_ZLIB):
		# type: (str, int) -> None
		if compression_type not in (COMPRESSION_TYPE_NONE, COMPRESSION_TYPE_ZLIB):
			raise ValueError('Unknown compression type {}'.format(compression_type))

		self._path = path
		self._data = open(path, "wb")
		self._dtype = None
		self._compression_type = compression_type
		self._current_data_block = []  # type: List[bytes]


	def append(self, value):
		# type: (Any) -> None
		"""
		Appends an array-like entry to the end of a column, converting the
		value to a 1-D array. Specific implementation depends on subclasses.

		Parameters:
			value (array-like): A NumPy ndarray or anything that can be cast to
				an array via np.asarray, including a scalar (i.e. 0D array).
		"""
		if self._data.closed:
			raise ValueError('I/O operation on closed file')


	def _write_block(self):
		# type: () -> None
		"""
		Compress and write the current block, if any. Specific implementation
		depends on subclasses.
		"""
		pass


	def _get_dtype_descr(self):
		# type: () -> bytes
		"""
		Get the description of the column data type in JSON format.
		"""
		assert self._dtype
		descr = self._dtype.descr
		if len(descr) == 1 and descr[0][0] == "":
			descr = descr[0][1]
		descr_json = json.dumps(descr, separators=(',', ':')).encode('utf-8')

		return descr_json


	def close(self):
		# type: () -> None
		"""
		Finish writing and close the column's data file. Idempotent.

		Notes
		-----
		Trying to append after closing will raise an error.
		"""
		if not self._data.closed:
			try:
				self._write_block()
				self._data.truncate()
			finally:
				self._data.close()


	def __del__(self):
		# type: () -> None
		"""
		Explicitly closes the output file once the instance is totally
		dereferenced.

		Notes
		-----
		This will lead to errors consequent of operating on a closed file if
		references to the output file (which ought to be private)
		exist. This is desirable, because such references should not be made.
		"""
		self.close()


class _FixedLengthColumn(_Column):
	"""
	Subclass for columns with fixed-length arrays for each row.

	The file generated by this class begins with a COLUMN_CHUNK_TYPE header
	chunk that specifies metadata (bytes per entry, entries per block, elements
	per entry, compression type, data type) for the column and is followed by
	one or more BLOCK_CHUNK_TYPE data chunks that contain the values of the
	actual rows. The data chunks are optionally compressed.
	"""
	def __init__(self, path, compression_type=COMPRESSION_TYPE_ZLIB):
		# type: (str, int) -> None
		super(_FixedLengthColumn, self).__init__(
			path, compression_type=compression_type)

		self._bytes_per_entry = 0
		self._entries_per_block = 0
		self._elements_per_entry = 0  # aka subcolumn count


	def append(self, value):
		# type: (Any) -> None
		"""
		Append a row to the fixed-length column.

		The first call to this method will define the column's NumPy dtype,
		element array size (subcolumns), and element size in bytes. Subsequent
		entries must be consistent in size.

		Parameters:
			value (array-like): A NumPy ndarray or anything that can be cast to
				an array via np.asarray, including a scalar (i.e. 0D array).
		"""
		super(_FixedLengthColumn, self).append(value)

		value = np.asarray(value, self._dtype)
		if value.nbytes == 0:
			# Replace an empty row with [0] to preserve num_rows without
			# special cases in the rest of the writer, reader, and downstream
			# code.
			value = np.zeros(1, value.dtype)
		data_bytes = value.tobytes()

		# First entry: Write the column header.
		if self._dtype is None:
			self._dtype = value.dtype

			self._bytes_per_entry = value.nbytes
			self._entries_per_block = int(np.ceil(BLOCK_BYTES_GOAL / value.nbytes))
			self._elements_per_entry = value.size
			descr_json = self._get_dtype_descr()

			chunk_header = CHUNK_HEADER.pack(
				COLUMN_CHUNK_TYPE,
				COLUMN_STRUCT.size + len(descr_json),
				)
			column_struct = COLUMN_STRUCT.pack(
				self._bytes_per_entry,     # I: bytes/entry, before packing, compression, and CHUNK_HEADER
				self._elements_per_entry,  # I: subcolumns
				self._entries_per_block,   # H: packed entries/block; the last block may be smaller
				self._compression_type,    # H: compression type code
				)

			self._data.write(chunk_header + column_struct + descr_json)

		# Later entry: Check size consistency.
		elif self._bytes_per_entry != len(data_bytes):
			raise VariableEntrySizeError(
				'Entry size in bytes, elements {} is inconsistent with {}'
				' for Table column {} which is not set up for variable lengths.'.format(
					(len(data_bytes), value.size),
					(self._bytes_per_entry, self._elements_per_entry),
					self._path))

		# Collect up an I/O block. Compress and write it when it's full.
		self._current_data_block.append(data_bytes)

		if len(self._current_data_block) >= self._entries_per_block:
			self._write_block()


	def _write_block(self):
		# type: () -> None
		"""
		Compress and write the current block, if any.
		"""
		if self._current_data_block:
			block_data = b''.join(self._current_data_block)

			if self._compression_type == COMPRESSION_TYPE_ZLIB:
				block_data = zlib.compress(block_data, ZLIB_LEVEL)

			block_header = CHUNK_HEADER.pack(BLOCK_CHUNK_TYPE, len(block_data))

			self._data.write(block_header + block_data)
			del self._current_data_block[:]


class _VariableLengthColumn(_Column):
	"""
	Subclass for columns with variable-length arrays for each row.

	The file generated by this class begins with a VARIABLE_COLUMN_CHUNK_TYPE
	header chunk that specifies metadata (compression type, data type) for the
	column and is followed by one or more pairs of ROW_SIZE_CHUNK_TYPE row size
	chunks that contain the number of array elements in each row and
	BLOCK_CHUNK_TYPE data chunks that contain the values from the actual rows.
	All BLOCK_CHUNK_TYPE chunks are preceded by a ROW_SIZE_CHUNK_TYPE chunk.
	The data chunks are optionally compressed.

	Since the TableReader uses np.nan as filler values for the values stored
	in variable-length columns, all values are forced to use the np.float64
	data type to avoid any unexpected behavior (np.nan is not defined for
	integer arrays).
	"""
	def __init__(self, path, compression_type=COMPRESSION_TYPE_ZLIB):
		# type: (str, int) -> None
		super(_VariableLengthColumn, self).__init__(path, compression_type)

		self._current_row_sizes_block = []  # type: List[int]
		self._remaining_bytes_in_block = BLOCK_BYTES_GOAL


	def append(self, value):
		# type: (Any) -> None
		"""
		Append a row to the variable-length column.

		The first call to this method will define the column's NumPy dtype.
		Subsequent entries can have different lengths, and the sizes of each
		entry are recorded separately.

		Parameters:
			value (array-like): A NumPy ndarray or anything that can be cast to
				an array via np.asarray, including a scalar (i.e. 0D array).
		"""
		super(_VariableLengthColumn, self).append(value)

		# All values must use the np.float64 data type
		value = np.asarray(value, VARIABLE_COLUMN_DATA_DTYPE)
		data_bytes = value.tobytes()

		# First entry: Write the column header.
		if self._dtype is None:
			self._dtype = value.dtype
			descr_json = self._get_dtype_descr()

			chunk_header = CHUNK_HEADER.pack(
				VARIABLE_COLUMN_CHUNK_TYPE,
				VARIABLE_COLUMN_STRUCT.size + len(descr_json)
				)
			column_struct = VARIABLE_COLUMN_STRUCT.pack(
				self._compression_type,  # H: compression type code
				)

			self._data.write(chunk_header + column_struct + descr_json)

		# Collect up an I/O block. Compress and write it when it's full.
		self._current_data_block.append(data_bytes)
		self._current_row_sizes_block.append(value.size)
		self._remaining_bytes_in_block -= value.nbytes

		if self._remaining_bytes_in_block <= 0:
			self._write_block()


	def _write_block(self):
		# type: () -> None
		"""
		Compress and write the current block, if any. All block data chunks
		are preceded by a row size chunk that specifies the number of array
		elements in each row.
		"""
		if self._current_data_block:
			row_size_data = np.array(
				self._current_row_sizes_block, ROW_SIZE_CHUNK_DTYPE).tobytes()
			row_size_header = CHUNK_HEADER.pack(
				ROW_SIZE_CHUNK_TYPE, len(row_size_data))

			block_data = b''.join(self._current_data_block)

			if self._compression_type == COMPRESSION_TYPE_ZLIB:
				block_data = zlib.compress(block_data, ZLIB_LEVEL)

			block_header = CHUNK_HEADER.pack(BLOCK_CHUNK_TYPE, len(block_data))

			self._data.write(row_size_header + row_size_data + block_header + block_data)
			self._remaining_bytes_in_block = BLOCK_BYTES_GOAL
			del self._current_row_sizes_block[:]
			del self._current_data_block[:]


class TableWriter(object):
	"""
	Generic live streaming output writer for NumPy ndarrays.

	NumPy can save and load arrays to disk.  This class provides a convenient
	interface to repeated appending of NumPy array data to an output file,
	which can be loaded one column at a time via the TableReader class.

	A Table has one or more named columns. Each (row x column) entry is a
	1-D NumPy array.

	Output file structure:

	<root directory> : Root path, provided by during instantiation.
		/attributes.json : A JSON file containing the attributes and metadata.
		/<column name> : A file per column containing the array data in an
				extensible sequence of chunks. See _Column.

	Parameters:
		path (str): Path to the directory to create.  All data will be saved
			within this directory.  It's OK if the directory already exists but
			not OK if it contains Table files since existing or concurrent
			Table files would confuse each other, so this will raise
			TableExistsError if its attributes.json file already exists.

	See also
	--------
	wholecell.io.tablereader.TableReader

	Notes
	-----
	The terms used in this class are adapted from the analogy of a spreadsheet
	or table.

	Each append() operation adds a new 'row' to the table containing one
	'entry' (or table cell) per column.  Each 'column' corresponds to a 'field'
	and holds a fixed data type for all of its entries. append() will convert
	scalar values and higher dimension arrays to 1D arrays. For fixed-length
	columns, each 1D array must be of fixed size, and each of the elements are
	called 'subcolumns'. For variable-length columns, the 1D arrays given for
	each row may have different sizes. Before appending any values to
	variable-length columns, the method set_variable_column_names() must first
	be called with the names of the columns as the argument.

	Both simple and structured ndarray dtypes are supported.  Structured arrays
	are a good way to work around the dimensional limitations.

	NumPy object arrays (e.g. arrays not of one or more standard dtypes) are
	not supported.  They might work, but under the hood will rely on pickling
	for saving and loading, which will be terribly slow.

	A Table also stores named 'attributes' meant for user-provided annotations
	that may be useful for downstream analysis or portability.  E.g. a list of
	element names for a vector-column.  Each attribute can be written only
	once -- attributes do not support time-series data.

	TODO (John): Test portability across machines (particularly, different
		operating systems).

	TODO (John): Consider writing all fields simultaneously
		as part of a structured array (i.e. a hybrid data type).
	"""

	def __init__(self, path):
		# type: (str) -> None
		self._path = filepath.makedirs(path)
		self._columns = None  # type: Optional[Dict[str, _Column]]
		self._variable_length_columns = set()  # type: Set[str]

		self._attributes = {}  # type: Dict[str, Any]
		self._attributes_filename = os.path.join(path, FILE_ATTRIBUTES)

		if (os.path.exists(self._attributes_filename)
				or os.path.exists(os.path.join(path, V2_DIR_COLUMNS))):
			raise TableExistsError('In {}'.format(self._path))

		# The column file's header chunk type mostly obviates the '_version'
		# attribute but writing the attributes file now lets the above check
		# also prevent competing TableWriters.
		self.writeAttributes(_version=VERSION)


	def append(self, **namesAndValues):
		# type: (**Any) -> None
		"""
		Write a new row of values, writing a 1-D NumPy array to each named
		column.

		The first call to this method will define the column names and dtypes.
		Subsequent calls will validate the names and types for consistency.

		Parameters:
			**namesAndValues: The array-like values to append to the ends
				of the named columns.

		Notes
		-----
		All fields must be provided every time this method is called.
		"""
		# First call - instantiate all columns
		if self._columns is None:
			self._columns = {
				name: _VariableLengthColumn(os.path.join(self._path, name))
				if name in self._variable_length_columns
				else _FixedLengthColumn(os.path.join(self._path, name))
				for name in namesAndValues
				}

		# Later calls - check for missing or unrecognized fields
		else:
			missingFields = six.viewkeys(self._columns) - six.viewkeys(namesAndValues)
			unrecognizedFields = six.viewkeys(namesAndValues) - six.viewkeys(self._columns)

			if missingFields:
				raise MissingFieldError(
					"Missing fields: {}".format(", ".join(missingFields))
					)

			if unrecognizedFields:
				raise UnrecognizedFieldError(
					"Unrecognized fields: {}".format(", ".join(unrecognizedFields))
					)

		for name, value in six.viewitems(namesAndValues):
			self._columns[name].append(value)


	def writeAttributes(self, **namesAndValues):
		# type: (**Any) -> None
		"""
		Writes JSON-serializable data.

		Oftentimes some additional data is needed to contextualize the data in
		the columns.  This method can be called to write JSON-serializable
		data (e.g. a list of strings) alongside the column data.

		Parameters:
			**namesAndValues: The named JSON-serializable attribute values.

		NOTE: TableWriter uses attribute names starting with "_" to hold its
		internal metadata.

		Notes
		-----
		This method can be called at any time, so long as an attribute name is
		not reused. That's because attributes are not designed for time-series
		data. [We could almost require all writeAttributes() to happen before
		append() except 'lengthSec' and 'endTime' attributes get written to
		'Main' at the end of each generation.]
		"""

		# check before modifying self._attributes
		sanitized = {}  # type: Dict[str, Any]

		for name, value in six.viewitems(namesAndValues):
			if name in self._attributes:
				raise AttributeAlreadyExistsError(
					"An attribute named '{}' already exists.".format(name)
					)

			try:
				if isinstance(value, np.ndarray):
					print("Warning - converting '{}' attribute from ndarray to list for JSON serialization.".format(name))
					value = value.tolist()

				json.dumps(value)  # test that it's JSON serializable

			except TypeError:
				raise AttributeTypeError(
					"Attribute '{}' value ({!r}) was not JSON serializable.".format(name, value)
					)

			sanitized[name] = value

		self._attributes.update(sanitized)
		filepath.write_json_file(self._attributes_filename, self._attributes, indent=1)


	def set_variable_length_columns(self, *column_names):
		# type: (*str) -> None
		"""
		Sets the names of columns that should have entries with variable
		lengths. This must be set before any values are appended to the column.

		Args:
			*column_names (str): Names of columns that should have entries with
				variable lengths.
		"""
		assert self._columns is None, 'Can set variable-length columns only before appending data rows'

		for name in column_names:
			self._variable_length_columns.add(name)


	def close(self):
		# type: () -> None
		"""
		Close the output files (columns).

		Notes
		-----
		Trying to append after closing will raise an error.

		"""
		if self._columns is not None:
			for column in six.viewvalues(self._columns):
				column.close()


	def __del__(self):
		# type: () -> None
		"""
		Close the output files once the instance is totally dereferenced.
		"""
		self.close()
