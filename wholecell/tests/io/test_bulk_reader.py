'''
Test cases for TableReader readColumn method performance

Requires:
	- sim_out_dir (str): first argument, must be a path to a simOut directory

Usage example:
	python wholecell/tests/io/test_bulk_reader.py out/sim_desc/wildtype_000000/000000/generation_000000/000000/simOut

TODO:
	- Add unit tests for TableReader and TableWriter

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/28/18
'''

from __future__ import absolute_import
from __future__ import division

import os
import sys
import time

import numpy as np

from wholecell.io.tablereader import TableReader

ITERS = 10
BLOCK_SIZE = 5000  # roughly number of proteins or RNA


def test_method(method, text):
	'''
	Tests a method for indexing into data from a reader

	Inputs:
		method (function): reads data and selects indices according to the
			desired method to be tested
		text (str): description of method
	'''

	start = time.time()
	for i in range(ITERS):
		counts = method()
	end = time.time()

	print('\t{}: {:.3f} s'.format(text, end - start))

	return counts

def test_old(reader, indices):
	'''
	Tests original readColumn method where all data is read and then indices
	are selected from the entire data matrix.

	Inputs:
		reader (TableReader object): file to read data from to test performance
		indices (numpy array of int): indices of data to select
	'''

	return test_method(lambda : reader.readColumn('counts')[:, indices].squeeze(), 'Old method')

def test_new_block(reader, indices):
	'''
	Tests new readColumn method where one chunk of data is read at each time point
	and then indices are selected.

	Inputs:
		reader (TableReader object): file to read data from to test performance
		indices (numpy array of int): indices of data to select
	'''

	return test_method(lambda : reader.readColumn('counts', indices), 'New method, block read')

def test_new_multiple(reader, indices):
	'''
	Tests new readColumn method where multiple reads are made at each time point
	for contiguous data.  Performance can be significantly slower than reading
	all of the data with a high number of seeks (dispersed indices).

	Inputs:
		reader (TableReader object): file to read data from to test performance
		indices (numpy array of int): indices of data to select
	'''

	return test_method(lambda : reader.readColumn('counts', indices, False), 'New method, multiple reads')

def test_functions(functions, text, reader, indices):
	'''
	Tests the given readColumn methods for performance and the same output

	Inputs:
		functions (list of functions): test functions from above to call
		text (str): description to display for the set of indices
		reader (TableReader object): file to read data from to test performance
		indices (numpy array of int): indices of data to select
	'''

	print('\n{} with {} iterations:'.format(text, ITERS))
	results = [f(reader, indices) for f in functions]
	for result in results[1:]:
		np.testing.assert_array_equal(results[0], result)

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

	# Sets of functions to test
	three_functions = [test_old, test_new_block, test_new_multiple]
	two_functions = [test_old, test_new_block]

	# Test reads
	## Single index
	indices = np.array([0])
	test_functions(three_functions, 'One index', bulk_molecules, indices)

	## First and last index
	indices = np.array([0, n_mols-1])
	test_functions(three_functions, 'First and last indices', bulk_molecules, indices)

	## Large block
	indices = np.array(range(BLOCK_SIZE))
	test_functions(three_functions, 'Block indices', bulk_molecules, indices)

	## 2 Large blocks
	indices = np.array(range(BLOCK_SIZE) + range(n_mols)[-BLOCK_SIZE:])
	test_functions(three_functions, 'Two blocks of indices', bulk_molecules, indices)

	## Dispersed reads - multiple reads method is slow so only test two methods
	indices = np.linspace(0, n_mols-1, BLOCK_SIZE, dtype=np.int64)
	test_functions(two_functions, 'Dispersed indices', bulk_molecules, indices)

	## Random reads - multiple reads method is slow so only test two methods
	indices = np.array(range(n_mols))
	np.random.shuffle(indices)
	indices = indices[:BLOCK_SIZE]
	test_functions(two_functions, 'Random indices', bulk_molecules, indices)

	## All indices
	indices = np.array(range(n_mols))
	test_functions(three_functions, 'All indices', bulk_molecules, indices)


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Need to supply sim out directory path as an argument")
		sys.exit(1)

	sim_out_dir = sys.argv[1]

	test_performance(sim_out_dir)
