
from __future__ import absolute_import, division, print_function

from chunk import Chunk
import os
import json
import numpy as np
import zlib

from wholecell.utils import filepath
from . import tablewriter as tw

__all__ = [
	"TableReader",
	]

SUPPORTED_COMPRESSION_TYPES = (tw.COMPRESSION_TYPE_NONE, tw.COMPRESSION_TYPE_ZLIB)


class TableReaderError(Exception):
	"""
	Base exception class for TableReader-associated exceptions.
	"""
	pass


class VersionError(TableReaderError):
	"""
	An error raised when the input files claim to be from a different format or
	version of the file specification.
	"""
	pass


class DoesNotExistError(TableReaderError):
	"""
	An error raised when a column or attribute does not seem to exist.
	"""
	pass


class _ColumnHeader(object):
	'''Column header info read from a Column file's first chunk.'''
	def __init__(self, chunk):
		if chunk.getname() != tw.COLUMN_CHUNK_TYPE:
			raise VersionError('Not a supported Column file format/version')

		header_struct = chunk.read(tw.COLUMN_STRUCT.size)
		(self.bytes_per_entry,
		self.elements_per_entry,
		self.entries_per_block,
		self.compression_type) = tw.COLUMN_STRUCT.unpack(header_struct)

		if self.compression_type not in SUPPORTED_COMPRESSION_TYPES:
			raise VersionError('Unsupported Column compression type {}'.format(
				self.compression_type))

		descr_json = chunk.read()
		descr = json.loads(descr_json)

		if isinstance(descr, basestring):
			self.dtype = str(descr)  # really the dtype.descr
		else:
			# numpy requires list-of-tuples-of-strings
			# TODO(jerry): Support triples?
			self.dtype = [(str(n), str(t)) for n, t in descr]


class TableReader(object):
	"""
	Reads output generated by TableWriter.

	Parameters:
		path (str): Path to the input location (a directory).

	See also
	--------
	wholecell.io.tablewriter.TableWriter
	wholecell.tests.io.measure_bulk_reader
	wholecell.tests.io.measure_zlib
	"docs/misc/byte_strings_to_2D_arrays.md"
	"""

	def __init__(self, path):
		self._path = path

		# Read the table's attributes file
		attributes_filename = os.path.join(path, tw.FILE_ATTRIBUTES)
		try:
			self._attributes = filepath.read_json_file(attributes_filename)

		except IOError as e:
			raise VersionError(
				"Could not read a table's attributes file ({})."
				" Version 2 tables are not supported."  # they could be...
				" Unzip all table files if needed.".format(attributes_filename), e)

		# Check if the table's version matches the expected version
		version = self._attributes['_version']
		if version != tw.VERSION:
			raise VersionError("Expected version {} but found version {}".format(
				tw.VERSION, version))

		# List the column file names. Ignore the 'attributes.json' file.
		self._columnNames = {p for p in os.listdir(path) if '.json' not in p}


	def readAttribute(self, name):
		"""
		Return an attribute value.

		Parameters:
			name (str): The attribute name.

		Returns:
			value (any): The attribute value, JSON-deserialized from a string.
		"""

		if name not in self._attributes:
			raise DoesNotExistError("No such attribute: {}".format(name))
		return self._attributes[name]


	def readColumn2D(self, name, indices=None):
		"""
		Load a full column (all rows). Each row entry is a 1-D NumPy array of
		subcolumns, so the result is a 2-D array row x subcolumn. This method
		can optionally read just a vertical slice of all those arrays -- the
		subcolumns at the given `indices`.

		The current approach collects up the compressed blocks, allocates the
		result array, then unpacks entries into it, keeping each decompressed
		block in memory only while copying into the result array. For a large
		column, this saves considerable RAM and a bit of time over collecting
		decompressed blocks then combining them via np.vstack().

		See docs/misc/byte_strings_to_2D_arrays.md for more design tradeoffs
		and their performance measurements.

		Parameters:
			name (str): The name of the column.
			indices (ndarray[int]): The subcolumn indices to select from each
				entry, or None to read in all data.

				If provided, this can give a performance boost for columns that
				are wide and tall.

				NOTE: The speed benefit might only be realized if the file is
				in the disk cache (i.e. the file has been recently read), which
				should typically be the case. This will still save RAM.

		Returns:
			ndarray: a writable 2-D NumPy array (row x subcolumn).

		TODO (jerry): Bring back the code to block-read `indices` of the data
			from uncompressed tables or after decompression, via seek + read or
			np.frombuffer(data, dtype, count, offset). It might be worthwhile
			only when header.entries_per_block == 1.

			The speed of various read methods is surprising and shape dependent.
			Techniques like `frombuffer(join(all_the_bytestrings))` or loop
			over `result[i, :] = frombuffer(byestring)` tended to take about as
			long in a simple test of BulkMolecules/counts but run slower in
			measure_bulk_reader.py. The differences are more pronounced for a
			smaller table like BulkMolecules/atpRequested.
		"""
		def decomp(raw_block):
			'''Decompress and unpack a raw block to an ndarray.'''
			data = decompressor(raw_block)
			entries = np.frombuffer(data, header.dtype).reshape(
				-1, header.elements_per_entry)
			if indices is not None:
				entries = entries[:, indices]
			return entries

		if name not in self._columnNames:
			raise DoesNotExistError("No such column: {}".format(name))

		entry_blocks = []

		# Read the header and read, decompress, and unpack all the blocks.
		with open(os.path.join(self._path, name), 'rb') as dataFile:
			chunk = Chunk(dataFile, align=False)
			header = _ColumnHeader(chunk)
			chunk.close()
			decompressor = (
				zlib.decompress if header.compression_type == tw.COMPRESSION_TYPE_ZLIB
				else lambda data_bytes: data_bytes)

			while True:
				try:
					chunk = Chunk(dataFile, align=False)
				except EOFError:
					break

				if chunk.getname() == tw.BLOCK_CHUNK_TYPE:
					raw = chunk.read()
					if len(raw) != chunk.getsize():
						raise EOFError('Data block cut short {}/{}'.format(
							len(raw), chunk.getsize()))
					entry_blocks.append(raw)

				chunk.close()  # skips to the next chunk

		# Decompress the last block to get its shape, then allocate the result.
		raw = None  # release the block ref
		last_entries = decomp(entry_blocks.pop())
		last_num_rows = last_entries.shape[0]
		num_rows = len(entry_blocks) * header.entries_per_block + last_num_rows
		num_subcolumns = header.elements_per_entry if indices is None else len(indices)
		result = np.zeros((num_rows, num_subcolumns), header.dtype)

		row = 0
		for raw in entry_blocks:
			entries = decomp(raw)
			additional_rows = entries.shape[0]
			result[row : (row + additional_rows)] = entries
			row += additional_rows

		result[row : (row + last_num_rows)] = last_entries
		return result


	def readColumn(self, name, indices=None):
		'''
		Read a column via readColumn2D() then squeeze() the resulting
		NumPy array into a 0D, 1D, or 2D array, depending on the number
		of rows and subcolumns it has.

		1 row x 1 subcolumn => 0D.

		n rows x 1 subcolumn or 1 row x m subcolumns => 1D.

		n rows x m subcolumns => 2D.

		Args:
			name (str): the column name.
			indices (ndarray[int]): The subcolumn indices to select from each
				entry, or None to read in all data. See readColumn2D().

		Returns:
			ndarray: A writable 0D, 1D, or 2D array.
		'''
		return self.readColumn2D(name, indices).squeeze()


	def allAttributeNames(self):
		"""
		Returns a list of all attribute names including Table metadata.
		"""
		return self._attributes.keys()


	def attributeNames(self):
		"""
		Returns a list of ordinary (client-provided) attribute names.
		"""
		names = [k for k in self._attributes.iterkeys() if not k.startswith('_')]
		return names


	def columnNames(self):
		"""
		Returns the names of all columns.
		"""
		return list(self._columnNames)


	def close(self):
		"""
		Does nothing.

		The TableReader keeps no files open, so this method does nothing.

		Notes
		-----
		TODO (John): Consider removing this method.  At the moment are usage is
			inconsistent, and gives the impression that it is actually
			beneficial or necessary.
		"""
		pass
