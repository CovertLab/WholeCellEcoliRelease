'''
Test cases for TableReader readColumn method performance

Requires:
	- sim_out_dir (str): first argument, must be a path to a simOut directory

Usage example:
	python wholecell/tests/io/measure_bulk_reader.py out/sim_desc/wildtype_000000/000000/generation_000000/000000/simOut

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/28/18
'''

from __future__ import absolute_import, division, print_function

import os
import sys
from typing import Callable, List

import numpy as np
from six.moves import range

from wholecell.io.tablereader import TableReader
from wholecell.utils.py3 import monotonic_seconds


ITERS = 10
BLOCK_SIZE = 5000  # roughly number of proteins or RNA


def test_method(method, text):
	# type: (Callable[[], np.ndarray], str) -> np.ndarray
	'''
	Tests a method for indexing into data from a reader

	Inputs:
		method (function): reads data and selects indices according to the
			desired method to be tested
		text (str): description of method
	'''

	counts = np.array([])
	start = monotonic_seconds()
	for i in range(ITERS):
		counts = method()
	end = monotonic_seconds()

	print('\t{}: {:.3f} s'.format(text, end - start))

	return counts

def test_old(reader, column, indices):
	'''
	Tests original readColumn method where all data is read and then subcolumn
	indices are selected from the entire data matrix.

	Inputs:
		reader (TableReader object): file to read data from to test performance
		column (str): the column name to read
		indices (numpy array of int): indices of data to select
	'''

	return test_method(
		lambda : reader.readColumn2D(column)[:, indices],
		'Old method, read all + select subcolumns')

def test_old_full(reader, column, indices):
	'''
	Tests original readColumn method where all data is read with no subcolumn
	selection by indices.

	Inputs:
		reader (TableReader object): file to read data from to test performance
		column (str): the column name to read
		indices (numpy array of int): indices of data to select [ignored]
	'''

	return test_method(
		lambda : reader.readColumn2D(column),
		'Old method, read all')

def test_new_block(reader, column, indices):
	'''
	Tests the new readColumn method that reads one I/O block per entry,
	then selects subcolumns by indices, then aggregates the rows.

	Inputs:
		reader (TableReader object): file to read data from to test performance
		column (str): the column name to read
		indices (numpy array of int): indices of data to select
	'''

	return test_method(
		lambda : reader.readColumn2D(column, indices),
		'New method, read subcolumns')

def test_functions(functions, text, reader, column, indices):
	# type: (List[Callable[[TableReader, str, np.ndarray], np.ndarray]], str, TableReader, str, np.ndarray) -> np.ndarray
	'''
	Tests the given readColumn methods for performance and the same output

	Inputs:
		functions (list of functions): test functions from above to call
		text (str): description to display for the set of indices
		reader (TableReader object): file to read data from to test performance
		indices (numpy array of int): indices of data to select
	'''

	print('\n{} "{}" with {} iterations:'.format(text, column, ITERS))
	results = [f(reader, column, indices) for f in functions]
	for result in results[1:]:
		np.testing.assert_array_equal(results[0], result)
	return results[1]

def test_performance(sim_out_dir):
	'''
	Performs tests on multiple index conditions to compare times of various methods.

	Inputs:
		sim_out_dir (str): directory of simulation output to read from
	'''

	# Bulk molecule information
	bulk_molecules = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	bulk_ids = bulk_molecules.readAttribute('objectNames')
	n_mols = len(bulk_ids)

	# Mass table
	mass = TableReader(os.path.join(sim_out_dir, 'Mass'))

	# Sets of functions to test
	two_functions = [test_old, test_new_block]
	three_functions = [test_old_full, test_old, test_new_block]

	# Test reads
	## Mass/time is a small column with just one int per entry
	indices = np.array([0])
	test_functions(three_functions, '1-element column', mass,
		'time', indices)

	## Single index
	indices = np.array([0])
	test_functions(two_functions, 'One index into', bulk_molecules,
		'counts', indices)

	## First and last index
	indices = np.array([0, n_mols-1])
	test_functions(two_functions, 'First and last indices into', bulk_molecules,
		'counts', indices)

	## Large block
	indices = np.array(list(range(BLOCK_SIZE)))
	test_functions(two_functions, 'Block indices into', bulk_molecules,
		'counts', indices)

	## 2 Large blocks
	indices = np.array(list(range(BLOCK_SIZE)) + list(range(n_mols))[-BLOCK_SIZE:])
	test_functions(two_functions, 'Two blocks of indices into', bulk_molecules,
		'counts', indices)

	## Dispersed reads
	indices = np.linspace(0, n_mols-1, BLOCK_SIZE, dtype=np.int64)
	test_functions(two_functions, 'Dispersed indices', bulk_molecules,
		'counts', indices)

	## Random reads
	indices = np.array(list(range(n_mols)))
	np.random.shuffle(indices)
	indices = indices[:BLOCK_SIZE]
	test_functions(two_functions, 'Random indices into', bulk_molecules,
		'counts', indices)

	## All indices, same large column as most of these tests
	indices = np.array(list(range(n_mols)))
	test_functions(three_functions, 'All indices into', bulk_molecules,
		'counts', indices)

	## All indices, narrow column
	# ASSUMES: The processNames attribute indicates how many subcolumns are in
	# the atpRequested table. Check after timing; beforehand could preload the
	# disk cache.
	n_processes = len(bulk_molecules.readAttribute('processNames'))
	indices = np.array(list(range(n_processes)))
	a = test_functions(three_functions, 'All indices, narrow column', bulk_molecules,
		'atpRequested', indices)
	if a.shape[1] != n_processes:
		raise ValueError('Expected {} "atpRequested" subcolumns; got {}'.format(
			n_processes, a.shape[1]))


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Need to supply sim out directory path as an argument")
		sys.exit(1)

	sim_out_dir = sys.argv[1]

	test_performance(sim_out_dir)
