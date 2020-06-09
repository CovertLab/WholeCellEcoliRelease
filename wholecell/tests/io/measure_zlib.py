'''
Space and time performance measurements of zlib on wcEcoli Table data:
	* zlib parameter variations,
	* compressing incrementally vs. separate blocks,
	* block sizes,
	* zlib flush modes, and
	* the presence of Z_SYNC_FLUSH 4-byte suffixes to skip and replace.

See https://www.bolet.org/~pornin/deflate-flush.html about flush modes.

Tests out how to do incremental zlib (de)compression. The key detail turns out
to be:
	* dc.decompress(dc.unconsumed_tail + bs)

Usage example:
	python wholecell/tests/io/measure_zlib.py \
	out/manual/wildtype_000000/000000/generation_000000/000000/simOut \
	BulkMolecules atpRequested

Also try it on:
	* BulkMolecules counts  # a big array
	* Mass time             # small, 8-byte elements

NOTE: These functions have no embedded blank lines so you can paste them into
an interactive Python session for experimentation.

Conclusions 2018-11-15:
	* Use compression level 6 (the default). Higher levels are very slow on
	  large blocks for little compression gain. Lower levels didn't save much
	  time.
	* It might be worth coming up with a heuristic that picks the compression
	  level for the data, e.g. based on the size of the first Table entry or
	  the compression time for the first packed block.
	* Packing small data items into a larger block improves compression speed
	  and compression ratio. 4K - 16K is good.
	* Incremental compression is only useful if you can't do packing. Flushing
	  the compressor costs output size and compression time but it's needed if
	  the receiver needs to decompress a message without waiting for later
	  messages. Negative WBITS (compression window size) saves 6 bytes.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2018-11-12
'''

from __future__ import absolute_import, division, print_function

from functools import reduce
import os
import sys
from time import clock
import zlib

import numpy as np
from six.moves import range

from wholecell.io.tablereader import TableReader


CHUNK_HEADER_SIZE = 8
FLUSH_SUFFIX = b'\x00\x00\xff\xff'


def timeit(f):
	'''Time the execution of f(). Return (seconds, f())'''
	c = clock
	start = c()
	result = f()
	end = c()
	return end - start, result

def compress_list_incremental(bytestrings, level=6, wbits=15, flush_each=False):
	'''zlib-compress a list of bytestrings at the given compression level to a
	list of compressed bytestrings.
	This flushes the compressor at the end, appending one more bytestring.
	If flush_each, this also does a Z_SYNC_FLUSH after each bytestring.
	Each output piece except for the last one, should end with
	FLUSH_SUFFIX which could be stripped off before transmission then
	restored before decompression.
	Q. Discard the empty strings in the output?
	'''
	co = zlib.compressobj(level, zlib.DEFLATED, wbits)
	pieces = []
	for bs in bytestrings:
		pieces.append(co.compress(bs))
		if flush_each:
			pieces.append(co.flush(zlib.Z_SYNC_FLUSH))
	pieces.append(co.flush(zlib.Z_FINISH))
	return pieces

def decompress_list_incremental(bytestrings, wbits=15):
	'''zlib-decompress a list of bytestrings.
	Q. Can there be empty strings in the output?
	'''
	dc = zlib.decompressobj(wbits)
	pieces = []
	for bs in bytestrings:
		pieces.append(dc.decompress(dc.unconsumed_tail + bs))
	pieces.append(dc.flush())
	return pieces

def compress_list_block(bytestrings, level=6):
	pieces = []
	for bs in bytestrings:
		pieces.append(zlib.compress(bs, level))
	return pieces

def decompress_list_block(bytestrings):
	pieces = []
	for bs in bytestrings:
		pieces.append(zlib.decompress(bs))
	return pieces

def ndarray_to_bytestrings(nd):
	'''Export a 2D NumPy array to a list of bytestrings.'''
	return [nd[i, :].tobytes() for i in range(nd.shape[0])]

def sum_len(bytestrings):
	'''Return the total length in a list of bytestrings.'''
	return reduce(lambda total, bs: total + len(bs), bytestrings, 0)

def join(bytestrings):
	'''Join a list of bytestrings.'''
	return b''.join(bytestrings)

def measure1_incremental(array, level, wbits, flushes):
	'''Test one incremental compression/decompression combo; return measurments.'''
	bytestrings = ndarray_to_bytestrings(array)
	input_size = sum_len(bytestrings)
	secs, compressed = timeit(
		lambda: compress_list_incremental(bytestrings, level, wbits, flushes))
	decompressed = decompress_list_incremental(compressed, wbits)
	assert join(bytestrings) == join(decompressed)
	count_pads = sum([b.endswith(FLUSH_SUFFIX) for b in compressed])
	count_chunks = count_pads + (not compressed[-1].endswith(FLUSH_SUFFIX))
	compressed_size = sum_len(compressed)
	compression_percent = compressed_size * 100 / input_size
	c8_size = compressed_size + CHUNK_HEADER_SIZE * count_chunks  # with chunk headers
	c8_percent = c8_size * 100 / input_size
	f = 'F' if flushes else '-'
	return level, wbits, f, input_size, compressed_size, secs, compression_percent, c8_percent, count_pads

def measure1_block(array, level):
	'''Test one block compression/decompression combo; return measurements.'''
	bytestrings = ndarray_to_bytestrings(array)
	input_size = sum_len(bytestrings)
	secs, compressed = timeit(
		lambda: compress_list_block(bytestrings, level))
	decompressed = decompress_list_block(compressed)
	assert join(bytestrings) == join(decompressed)
	compressed_size = sum_len(compressed)
	compression_percent = compressed_size * 100 / input_size
	c8_size = compressed_size + CHUNK_HEADER_SIZE * len(compressed)  # with chunk headers
	c8_percent = c8_size * 100 / input_size
	return level, input_size, compressed_size, secs, compression_percent, c8_percent

def measure_incremental(array, table_name='', column_name=''):
	'''Measure incremental compression on the given numpy array w/various
	parameters. This incrementally compresses one row at a time, repeating the
	computation at multiple compression levels, with multiple WBITS parameters,
	and with and without flushing the compressor after each row.

	Parameters:
		array (ndarray): the numpy array to compress
		table_name (str): the source table name to label the printout
		column_name (str): the source column name to label the printout
	'''
	print('measure_incremental {}/{}'.format(table_name, column_name))
	print('{:<4} {:3} {:3} {:1} {:>9} {:>8} {:>8} {:>7} {:>7} {:>4}'.format(
		'FMT', 'LVL', 'WBT', 'F', 'I BYTES', 'C BYTES', 'SECS', 'COMP %', 'C8 %', 'PADS'))
	for level in (6, 9):
		for wbits in (15, -15):
			for flush in (False, True):
				tup = measure1_incremental(array, level, wbits, flush)
				print('zlib {:<3d} {:3d} {:1} {:9d} {:8d} {:8.4f} {:6.1f}% {:6.1f}% {:4d}'.format(*tup))

def measure_block(array):
	'''Measure block compression on the given numpy array w/various parameters.'''
	print('{:<4} {:3} {:>9} {:>9} {:>8} {:>8} {:>7} {:>7}'.format(
		'FMT', 'LVL', 'ROW BYTES', 'I BYTES', 'C BYTES', 'SECS', 'COMP %', 'C8 %'))
	row_bytes = array[0, :].nbytes
	levels = (4, 5, 6, 7, 8, 9) if array[0].nbytes < 400000 else (5, 6)
	for level in levels:
		tup = measure1_block(array, level)
		tup = (tup[0], row_bytes) + tup[1:]
		print('zlib {:<3d} {:9d} {:9d} {:8d} {:8.4f} {:6.1f}% {:6.1f}%'.format(*tup))

def measure_tables(basepath, table_columns):
	for table, columns in table_columns:
		reader = TableReader(os.path.join(basepath, table))
		for column in columns:
			array = reader.readColumn2D(column)
			# TODO(jerry): ...

def time_decompressor(bytestrings, level):
	'''Measure the block decompression time at different compression levels.'''
	# Results: Decompression time is about the same at levels 4 - 9, and maybe
	# 30% faster only with large data at levels 0 - 3.
	pieces = compress_list_block(bytestrings, level)
	results = timeit(lambda: decompress_list_block(pieces))
	return results[0]

def measure_block_subranges(array, table_name='', column_name=''):
	'''Measure block compression on aggregated subranges of the given array.
	This is to help with choosing compression parameters.

	Parameters:
		array (ndarray): the numpy array to compress
		table_name (str): the source table name to label the printout
		column_name (str): the source column name to label the printout
	'''
	print('measure_block_subranges {}/{}'.format(table_name, column_name))
	rows, columns = array.shape
	rows_log2 = int(np.log2(rows))
	max_ = 2 ** rows_log2
	power_of_2_array = array[:max_, :]  # limited to 2**max rows for this test
	for power in range(0, rows_log2 + 1):
		rows_per_block = 2 ** power
		data = power_of_2_array.reshape(max_ // rows_per_block, -1)
		measure_block(data)
	return power_of_2_array

def measure_performance_metrics(sim_out_dir, table_name, column_name):
	reader = TableReader(os.path.join(sim_out_dir, table_name))
	array = reader.readColumn2D(column_name)
	print('')
	array2 = measure_block_subranges(array, table_name, column_name)
	print('')
	measure_incremental(array2, table_name, column_name)
	print('''
FMT = compression format
LVL = compression level
ROW BYTES = bytes/entry (a column's row)
I BYTES = total number of Input bytes processed
C BYTES = total number of Compressed bytes output
SECS = compression time
COMP % = compression ratio in %
C8 % = compression ratio in % including 8-byte chunk headers
PADS = number of 'x00x00xffxff' zlib padding suffixes
	(The pad suffixes don't need to be stored if they're removed after
	compression and reinserted before decompression.)
	''')

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print("Args needed: simOut_directory, table_name, column_name")
		sys.exit(1)

	measure_performance_metrics(*sys.argv[1:4])
