"""Utilities for working with tables"""


from __future__ import absolute_import, division, print_function
from os import path

import numpy as np

from wholecell.io.tablereader import TableReader


SUBCOLS_NAME_MAP_ATTR = "subcolumns"


def read_subcolumn(sim_out_dir, table, column, subcolumn_name):
	# type: (str, str, str, str) -> np.ndarray
	"""Read in a subcolumn from a table by name

	Each column of a table is a 2D matrix. The SUBCOLS_NAMES_ATTR attribute
	defines a map from column name to a name for an attribute that
	stores a list of names such that the i-th name describes the i-th
	subcolumn.

	Arguments:
		sim_out_dir: Path to the simulation output directory.
		table: Name of the table.
		column: Name of the column.
		subcolumn_name: Name of the ID or object associated with the
			desired subcolumn.

	Returns:
		The subcolumn, as a 1-dimensional array.
	"""
	reader = TableReader(path.join(sim_out_dir, table))
	subcol_name_map = reader.readAttribute(SUBCOLS_NAME_MAP_ATTR)
	subcols = reader.readAttribute(subcol_name_map[column])
	index = subcols.index(subcolumn_name)
	return reader.readColumn2D(column, [index])[:, 0]
