#!/usr/bin/env python2


"""Tests for tablereader_utils utility functions"""


from __future__ import absolute_import, division, print_function
from os import path
from shutil import rmtree
from tempfile import mkdtemp
from typing import Dict, List
import unittest

import numpy as np

from wholecell.io.tablewriter import TableWriter
from prototypes.behavior_metrics.tablereader_utils import read_subcolumn


TABLE_NAME = "table"
COLUMN_NAME = "col"
COLUMN_2_NAME = "col2"
SUBCOL_NAMES = ["A", "B", "C", "D", "E", "F"]
SUBCOL_2_NAMES = ["Z", "Y", "X", "W", "V", "U"]


class TestReadSubcolumn(unittest.TestCase):

	def setUp(self):
		self.temp_dir = mkdtemp()
		self.tablewriter = TableWriter(
			path.join(self.temp_dir, TABLE_NAME)
		)

	def tearDown(self):
		rmtree(self.temp_dir)

	def test_read_normal_subcolumn(self):
		vals = np.random.random((5, 6))
		for val in vals:
			self.tablewriter.append(**{COLUMN_NAME: val})
		self._write_subcolumn_metadata({COLUMN_NAME: SUBCOL_NAMES})
		self.tablewriter.close()

		for i in range(vals.shape[1]):
			subcol = read_subcolumn(
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

		subcol = read_subcolumn(
			self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[0]
		)
		np.testing.assert_array_equal(vals, subcol)

	def test_read_sole_row(self):
		val = np.random.random()
		self.tablewriter.append(**{COLUMN_NAME: val})
		self._write_subcolumn_metadata({COLUMN_NAME: SUBCOL_NAMES[0:1]})
		self.tablewriter.close()

		subcol = read_subcolumn(
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
			subcol = read_subcolumn(
				self.temp_dir, TABLE_NAME, COLUMN_NAME, SUBCOL_NAMES[i]
			)
			np.testing.assert_array_equal(
				vals[:, 0, i], subcol
			)
			subcol = read_subcolumn(
				self.temp_dir, TABLE_NAME, COLUMN_2_NAME, SUBCOL_2_NAMES[i]
			)
			np.testing.assert_array_equal(
				vals[:, 1, i], subcol
			)

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
			col_name: str(i) for i, col_name in
			enumerate(subcolumns_labels.keys())
		}
		attributes = {
			str(i): lst for i, lst in
			enumerate(subcolumns_labels.values())
		}
		attributes["subcolumns"] = subcolumns_key
		self.tablewriter.writeAttributes(**attributes)


if __name__ == "__main__":
	unittest.main()
