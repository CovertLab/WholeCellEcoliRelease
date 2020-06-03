from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile
from typing import List, Dict, Union
import unittest

import numpy as np
import numpy.testing as npt


# def noop_decorator(fcn):
# 	return fcn
#
# # Workaround for @profile only defined when running under kernprof:
# __builtins__.setdefault('profile', noop_decorator)

from wholecell.utils import filepath
from wholecell.io.tablereader import (TableReader, DoesNotExistError,
	VariableLengthColumnError)
from wholecell.io.tablewriter import (BLOCK_BYTES_GOAL,
	TableWriter, MissingFieldError, TableExistsError, UnrecognizedFieldError,
	VariableEntrySizeError, AttributeAlreadyExistsError, AttributeTypeError,
	V2_DIR_COLUMNS)


COLUMNS = 'x y z theta'.split()
DATA = {key: np.arange(10.0) + ord(key[0]) for key in COLUMNS}

# For subcolumns read test
TABLE_NAME = "table"
COLUMN_NAME = "col"
COLUMN_2_NAME = "col2"
VARIABLE_LENGTH_COLUMN_NAME = "vcol"
SUBCOL_NAMES = ["A", "B", "C", "D", "E", "F"]
SUBCOL_2_NAMES = ["Z", "Y", "X", "W", "V", "U"]


# TODO(jerry): Test structured dtypes.

class Test_TableReader_Writer(unittest.TestCase):
	def setUp(self):
		self.test_dir = None
		self.table_path = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		'''Create a temp test dir. TableWriter must make its own subdir.'''
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()
			self.table_path = os.path.join(self.test_dir, 'Main')

	def test_basic(self):
		'''Test a table with some float arrays and no attributes.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(**DATA)

		d2 = {key: -10 * value for key, value in DATA.iteritems()}
		writer.append(**d2)

		no_theta = dict(d2)
		del no_theta['theta']
		with self.assertRaises(MissingFieldError):
			writer.append(**no_theta)

		with self.assertRaises(UnrecognizedFieldError):
			writer.append(JUNK=np.arange(4), **d2)

		writer.close()

		with self.assertRaises(ValueError):
			writer.append(**DATA)  # writer is closed

		# --- Read ---
		reader = TableReader(self.table_path)
		self.assertEqual([], reader.attributeNames())
		self.assertEqual({'_version'}, set(reader.allAttributeNames()))
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		with self.assertRaises(DoesNotExistError):
			reader.readAttribute('x')

		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)  # the basic readColumn() case
		expected = np.vstack((DATA[column_name], d2[column_name]))
		npt.assert_array_equal(expected, actual)
		self.assertEqual(2, actual.ndim)

		with self.assertRaises(DoesNotExistError):
			reader.readColumn('JUNK')

		reader.close()

	def test_readColumn2D(self):
		'''Test readColumn2D() vs. readColumn() array squeezing.
		readColumn2D() should always return a 2D array.
		readColumn() may return a 0D, 1D, or 2D array, depending on the data.
		Cases: {1, 2} rows x {scalar, 1-element, 1D, 2D} x {None, 1, 2} indices.
		TODO(jerry): Test 0 rows?
		'''
		self.make_test_dir()

		data = {
			'scalar': 1,
			'1_element': np.array([100]),
			'1D': np.arange(4, dtype=np.int8),
			'2D': np.arange(12).reshape(4, 3),
			}
		index1 = np.array([0])
		index2 = np.array([0, 1])

		# === Write 1 row ===
		writer = TableWriter(self.table_path)
		writer.append(**data)
		writer.close()

		# --- Read 1 row with readColumn2D() ---
		reader = TableReader(self.table_path)
		self.assertEqual(2, reader.readColumn2D('scalar').ndim)
		self.assertEqual(2, reader.readColumn2D('1_element').ndim)
		self.assertEqual(2, reader.readColumn2D('1D').ndim)
		self.assertEqual(2, reader.readColumn2D('2D').ndim)

		self.assertEqual(2, reader.readColumn2D('scalar', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('1_element', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('1D', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('2D', index1).ndim)

		with self.assertRaises(IndexError):
			self.assertEqual(2, reader.readColumn2D('scalar', index2).ndim)
		with self.assertRaises(IndexError):
			self.assertEqual(2, reader.readColumn2D('1_element', index2).ndim)
		self.assertEqual(2, reader.readColumn2D('1D', index2).ndim)
		self.assertEqual(2, reader.readColumn2D('2D', index2).ndim)

		# --- Read 1 row with readColumn() ---
		self.assertEqual(0, reader.readColumn('scalar').ndim)
		self.assertEqual(0, reader.readColumn('1_element').ndim)
		self.assertEqual(1, reader.readColumn('1D').ndim)
		self.assertEqual(1, reader.readColumn('2D').ndim)

		self.assertEqual(0, reader.readColumn('scalar', index1).ndim)
		self.assertEqual(0, reader.readColumn('1_element', index1).ndim)
		self.assertEqual(0, reader.readColumn('1D', index1).ndim)
		self.assertEqual(0, reader.readColumn('2D', index1).ndim)

		self.assertEqual(1, reader.readColumn('1D', index2).ndim)
		self.assertEqual(1, reader.readColumn('2D', index2).ndim)

		# Test that the reader returns writeable arrays.
		reader.readColumn2D('2D')[0, 0] = 0
		reader.readColumn2D('2D', index2)[0, 0] = 0

		# === Write 2 rows ===
		table2_path = os.path.join(self.test_dir, 'Segundo')
		writer = TableWriter(table2_path)
		writer.append(**data)
		writer.append(**data)
		writer.close()

		# --- Read 2 rows with readColumn2D() ---
		reader = TableReader(table2_path)
		self.assertEqual(2, reader.readColumn2D('scalar').ndim)
		self.assertEqual(2, reader.readColumn2D('1_element').ndim)
		self.assertEqual(2, reader.readColumn2D('1D').ndim)
		self.assertEqual(2, reader.readColumn2D('2D').ndim)

		self.assertEqual(2, reader.readColumn2D('scalar', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('1_element', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('1D', index1).ndim)
		self.assertEqual(2, reader.readColumn2D('2D', index1).ndim)

		self.assertEqual(2, reader.readColumn2D('1D', index2).ndim)
		self.assertEqual(2, reader.readColumn2D('2D', index2).ndim)

		# --- Read 2 rows with readColumn() ---
		self.assertEqual(1, reader.readColumn('scalar').ndim)
		self.assertEqual(1, reader.readColumn('1_element').ndim)
		self.assertEqual(2, reader.readColumn('1D').ndim)
		self.assertEqual(2, reader.readColumn('2D').ndim)

		self.assertEqual(1, reader.readColumn('scalar', index1).ndim)
		self.assertEqual(1, reader.readColumn('1_element', index1).ndim)
		self.assertEqual(1, reader.readColumn('1D', index1).ndim)
		self.assertEqual(1, reader.readColumn('2D', index1).ndim)

		self.assertEqual(2, reader.readColumn('1D', index2).ndim)
		self.assertEqual(2, reader.readColumn('2D', index2).ndim)

	def test_attributes(self):
		'''Test a table with attributes and no columns.'''
		def check_attributes(attribute_dict):
			for k, v in attribute_dict.iteritems():
				self.assertEqual(attribute_dict[k], reader.readAttribute(k))

		self.make_test_dir()

		d1 = dict(mercury=1, venus=2, earth=3, mars=4)
		d2 = dict(jupiter=[50, 60, 70], saturn='Saturn')
		d3 = dict(uranus=700.0, neptune=800.5)
		keys = set(d1.keys() + d2.keys() + d3.keys())

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.writeAttributes(**d1)
		writer.writeAttributes(**d2)
		writer.writeAttributes(**d3)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		self.assertEqual([], reader.columnNames())
		self.assertEqual(keys, set(reader.attributeNames()))
		self.assertEqual(len(keys), len(reader.attributeNames()))

		check_attributes(d1)
		check_attributes(d3)
		check_attributes(d2)
		check_attributes(d1)
		check_attributes(d2)

		with self.assertRaises(DoesNotExistError):
			reader.readColumn('JUNK')

		reader.close()

	def test_array_dimensions(self):
		'''Test 0D, 2D, and non-uniform array lengths, also automatic int to
		float value conversion.
		'''
		self.make_test_dir()
		d0 = {key: 19 for key in DATA.iterkeys()}
		d2 = {key: value.reshape(2, -1) for key, value in DATA.iteritems()}
		d3 = {key: value[2:] for key, value in DATA.iteritems()}

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(**DATA)  # 1-D float arrays in table row 0

		with self.assertRaises(VariableEntrySizeError):
			writer.append(**d0)  # 0-D arrays: inconsistent entry sizes

		writer.append(**d2)  # 2-D arrays in table row 2

		with self.assertRaises(VariableEntrySizeError):
			writer.append(**d3)  # narrower 1-D arrays than row 0

		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)
		self.assertEqual(1, actual[1].ndim)  # 2-D d2 values reshaped to 1-D
		expected = np.vstack((DATA[column_name], DATA[column_name]))
		npt.assert_array_equal(expected, actual)

	def test_scalars(self):
		'''Test the case where all rows contain scalars, where readColumn()
		returns a 1-D array instead of a 2-D array!
		'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(x=20)  # scalar int
		writer.append(x=21)
		writer.append(x=22)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual = reader.readColumn('x')
		self.assertEqual(1, actual.ndim)
		self.assertEqual((3,), actual.shape)
		npt.assert_array_equal([20, 21, 22], actual)

	def test_1_element_arrays(self):
		'''Test the case where all rows contain 1-element arrays, where
		readColumn() returns a 1-D array instead of a 2-D array!
		'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(x=[20])
		writer.append(x=[21])
		writer.append(x=[22])
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual = reader.readColumn('x')
		self.assertEqual(1, actual.ndim)
		self.assertEqual((3,), actual.shape)
		npt.assert_array_equal([20, 21, 22], actual)

	def test_0_element_arrays(self):
		'''Test the case where all rows contain 0-element arrays.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(val=[])  # an empty row becomes [0] to preserve num_rows
		writer.append(val=[])
		writer.append(val=[])
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual = reader.readColumn2D('val')
		self.assertEqual(2, actual.ndim)
		self.assertEqual((3, 1), actual.shape)
		npt.assert_array_equal(np.zeros((3, 1)), actual)

	def test_big_entries(self):
		'''Test entries that are larger than BLOCK_BYTES_GOAL.'''
		self.make_test_dir()
		ints = np.arange(9001, dtype=np.int16)
		floats = np.arange(9002, dtype=np.float16)
		d0 = {'ManyInts': ints, 'ManyFloats': floats}

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(**d0)
		writer.append(**d0)
		writer.append(**d0)
		writer.append(**d0)
		writer.append(**d0)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual_ints = reader.readColumn2D('ManyInts')
		actual_floats = reader.readColumn2D('ManyFloats')
		npt.assert_array_equal(np.vstack(5 * [ints]), actual_ints)
		npt.assert_array_equal(np.vstack(5 * [floats]), actual_floats)

	def test_many_entries(self):
		'''Test enough entries to pack into many blocks.'''
		self.make_test_dir()
		ints = np.arange(500, dtype=np.int8)
		floats = np.arange(500, dtype=np.float16)
		d0 = {'ManyInts': ints, 'ManyFloats': floats}

		rows = (BLOCK_BYTES_GOAL + ints.nbytes - 1) // ints.nbytes * 11 + 2

		# --- Write ---
		writer = TableWriter(self.table_path)
		for _ in xrange(rows):
			writer.append(**d0)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual_ints = reader.readColumn2D('ManyInts')
		actual_floats = reader.readColumn2D('ManyFloats')
		npt.assert_array_equal(np.vstack(rows * [ints]), actual_ints)
		npt.assert_array_equal(np.vstack(rows * [floats]), actual_floats)

	def test_variable_length_columns(self):
		'''Test variable-length columns.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.set_variable_length_columns(VARIABLE_LENGTH_COLUMN_NAME)
		writer.append(**{VARIABLE_LENGTH_COLUMN_NAME: np.arange(1, dtype=np.int32)})
		writer.append(**{VARIABLE_LENGTH_COLUMN_NAME: np.arange(2)})
		writer.append(**{VARIABLE_LENGTH_COLUMN_NAME: []})
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		variable_column = reader.readColumn(VARIABLE_LENGTH_COLUMN_NAME)

		npt.assert_array_equal(
			np.array([[0, np.nan], [0, 1], [np.nan, np.nan]]),
			variable_column)

		# Subcolumns cannot be accessed for variable-length columns
		with self.assertRaises(VariableLengthColumnError):
			reader.readColumn(VARIABLE_LENGTH_COLUMN_NAME, indices=1)

	def test_long_variable_length_columns(self):
		'''Test long variable-length columns that span multiple blocks.'''
		self.make_test_dir()

		rows = []
		row_lengths = []
		while len(rows) == 0 or np.concatenate(rows).nbytes < 10*BLOCK_BYTES_GOAL:
			row_length = np.random.randint(5, 100)
			rows.append(np.arange(row_length))
			row_lengths.append(row_length)

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.set_variable_length_columns(VARIABLE_LENGTH_COLUMN_NAME)
		for row in rows:
			writer.append(**{VARIABLE_LENGTH_COLUMN_NAME: row})
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		variable_column = reader.readColumn(VARIABLE_LENGTH_COLUMN_NAME)

		actual_column = np.full((len(rows), np.array(row_lengths).max()), np.nan)
		for i, row in enumerate(rows):
			actual_column[i, :row_lengths[i]] = row

		npt.assert_array_equal(actual_column, variable_column)

	def test_empty_variable_length_columns(self):
		'''Test variable-length columns with empty rows.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.set_variable_length_columns(VARIABLE_LENGTH_COLUMN_NAME)
		for i in range(1000):
			writer.append(
				**{VARIABLE_LENGTH_COLUMN_NAME: np.array([])})
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		variable_column = reader.readColumn(VARIABLE_LENGTH_COLUMN_NAME)

		npt.assert_array_equal(np.zeros((1000, 0)), variable_column)

	def test_path_clash(self):
		'''Test two TableWriters trying to write to the same directory.'''
		self.make_test_dir()

		TableWriter(self.table_path)

		with self.assertRaises(TableExistsError):
			TableWriter(self.table_path)

	def test_v2_clash(self):
		'''Test a TableWriter trying to write to a V2 table's directory.'''
		self.make_test_dir()
		filepath.makedirs(self.table_path, V2_DIR_COLUMNS)

		with self.assertRaises(TableExistsError):
			TableWriter(self.table_path)

	def test_attribute_errors(self):
		'''Test errors adding attributes.'''
		self.make_test_dir()

		writer = TableWriter(self.table_path)

		with self.assertRaises(AttributeAlreadyExistsError):
			writer.writeAttributes(_version=-1)
		writer.writeAttributes(heading=90.0)  # attributes should still work

		with self.assertRaises(AttributeAlreadyExistsError):
			writer.writeAttributes(heading=-90.0)
		writer.writeAttributes(rate=60.0, stars='*****')

		with self.assertRaises(AttributeTypeError):
			writer.writeAttributes(object=self)
		writer.writeAttributes(speed=186000.0, list=['a', 'b', 'c'])

		with self.assertRaises(AttributeTypeError):
			writer.writeAttributes(fcn=lambda: 1)
		writer.writeAttributes(x=1, y=2, z=3, red='green')

		with self.assertRaises(AttributeTypeError):
			writer.writeAttributes(fcn=lambda: 1, object=self, alpha='zed')
		writer.writeAttributes(alpha='ZED')  # not previously written




class TestReadSubcolumn(unittest.TestCase):

	def setUp(self):
		self.temp_dir = tempfile.mkdtemp()
		self.tablewriter = TableWriter(
			os.path.join(self.temp_dir, TABLE_NAME)
		)

	def tearDown(self):
		shutil.rmtree(self.temp_dir)

	def test_read_normal_subcolumn(self):
		vals = np.random.random((5, 6))
		for val in vals:
			self.tablewriter.append(**{COLUMN_NAME: val})
		self._write_subcolumn_metadata({COLUMN_NAME: SUBCOL_NAMES})
		self.tablewriter.close()

		for i in range(vals.shape[1]):
			subcol = self._read_subcolumn(
				self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[i]
			)
			np.testing.assert_array_equal(
				vals[:, i], subcol
			)

	def test_read_sole_subcolumn(self):
		vals = np.random.random(5)
		for val in vals:
			self.tablewriter.append(**{COLUMN_NAME: val})
		self._write_subcolumn_metadata({COLUMN_NAME: SUBCOL_NAMES[0:1]})
		self.tablewriter.close()

		subcol = self._read_subcolumn(
			self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[0]
		)
		np.testing.assert_array_equal(vals, subcol)

	def test_read_sole_row(self):
		val = np.random.random()
		self.tablewriter.append(**{COLUMN_NAME: val})
		self._write_subcolumn_metadata({COLUMN_NAME: SUBCOL_NAMES[0:1]})
		self.tablewriter.close()

		subcol = self._read_subcolumn(
			self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[0]
		)
		np.testing.assert_array_equal(np.array(val), subcol)

	def test_read_multiple_columns(self):
		vals = np.random.random((5, 2, 6))
		for val in vals:
			self.tablewriter.append(
				**{COLUMN_NAME: val[0], COLUMN_2_NAME: val[1]}
			)
		self._write_subcolumn_metadata(
			{COLUMN_NAME: SUBCOL_NAMES, COLUMN_2_NAME: SUBCOL_2_NAMES}
		)
		self.tablewriter.close()

		for i in range(vals.shape[2]):
			subcol = self._read_subcolumn(
				self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[i]
			)
			np.testing.assert_array_equal(
				vals[:, 0, i], subcol
			)
			subcol = self._read_subcolumn(
				self.temp_dir, TABLE_NAME, COLUMN_2_NAME, SUBCOL_2_NAMES[i]
			)
			np.testing.assert_array_equal(
				vals[:, 1, i], subcol
			)

	def _read_subcolumn(self, sim_out, table, column, subcolumn):
		reader = TableReader(os.path.join(sim_out, table))
		return reader.readSubcolumn(column, subcolumn)

	def _write_subcolumn_metadata(self, subcolumns_labels):
		# type: (Dict[str, List[str]]) -> None
		"""Record subcolumn metadata in table attributes

		Precondition:
			self.tablewriter must be ready to write
		Postcondition:
			The table will have an attribute `subcolumns` that stores a
			dictionary mapping from column names to some arbitrary name.
			For each such arbitrary name, the table will have another
			attribute of that name which stores the list of the names
			for the subcolumns of the associated column, in the order in
			which the subcolumns are indexed.
		Arguments:
			subcolumns_labels: A dictionary where for each key-value
				pair, the key is the name of a column in the table and
				the value is a list of the labels for each of that
				column's subcolumns such that the i-th subcolumn is
				named by the i-th element of the list.
		"""
		subcolumns_key = {
			col_name: str(i)
			for i, col_name in enumerate(subcolumns_labels.keys())
			}  # type: Dict[str, str]
		attributes = {
			str(i): lst
			for i, lst in enumerate(subcolumns_labels.values())
			}  # type: Dict[str, Union[List[str], Dict[str, str]]]
		attributes["subcolumns"] = subcolumns_key
		self.tablewriter.writeAttributes(**attributes)

