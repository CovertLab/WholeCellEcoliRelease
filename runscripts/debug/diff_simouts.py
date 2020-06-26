'''
Compare whether two simulation runs produced identical output and, if not,
print info on their differences. This also returns a shell exit status code.

This aim is to check that a code change such as a simulation speedup didn't
accidentally change the output. This does not support a tolerance between nor
check for semantic equivalence.

Example command line:

    diff_simouts.py out/manual/wildtype_000000/000000/generation_000000/000000/simOut \
    	out/experiment/wildtype_000000/000000/generation_000000/000000/simOut

See the `-h` command line usage help.
'''

from __future__ import absolute_import, division, print_function

import numpy as np
import os
import os.path as op
from pprint import pprint
import re
import sys

from wholecell.io.tablereader import TableReader, VersionError


WHITESPACE = re.compile(r'\s+')

class Repr(object):
	'''A Repr prints the given repr_ string without quotes and is not equal
	to any other value.
	'''
	def __init__(self, repr_):
		self.repr_ = repr_

	def __repr__(self):
		return self.repr_

def is_repr(value):
	'''Returns True if the value is a Repr object.'''
	return isinstance(value, Repr)

def elide(value, max_=100):
	'''Return a value with the same repr but elided if it'd be longer than max.'''
	repr_ = repr(value)
	if len(repr_) > max_:
		return Repr(repr_[:max_] + '...')
	return value

def compare_arrays(array1, array2):
	'''Compare two ndarrays, checking the shape and all elements, allowing for
	NaN values and non-numeric values.

	Args:
		array1 (np.ndarray): array to compare
		array2 (np.ndarray): array to compare
	Return:
		description (Repr): a "representation" of the array mismatch, or '' if
			the arrays match
	'''
	try:
		np.testing.assert_array_equal(array1, array2)
		return ''
	except AssertionError as e:
		return elide(Repr(WHITESPACE.sub(' ', e.args[0]).strip()))

def open_table(simout_dir, subdir):
	'''Try to open a Table in the given simOut/subdir.

	Args:
		simout_dir (str): a simOut directory
		subdir (str): a subdirectory name within simout_dir
	Returns:
		tuple (tuple[TableReader, Repr): a tuple containing a TableReader or
			None if unsuccessful, and a descriptive Repr
	'''
	table_path = op.join(simout_dir, subdir)

	try:
		table_reader = TableReader(table_path)
		table_message = Repr('<valid Table>')
	except VersionError:
		# Note: This could use the exception message but it's verbose.
		table_reader = None
		table_message = Repr('<not a Table>' if op.isdir(table_path) else '<absent subdir>')

	return table_reader, table_message

def diff_subdirs(subdir, simout_dir1, simout_dir2):
	'''Use the Table contents of the named subdir in `simout_dir1` as a
	reference point to diff the corresponding Table contents in `simout_dir2`.

	Args:
		subdir (str): the name of a subdir, hopefully containing a Table in
			both simOut dirs
		simout_dir1 (str): the reference simOut directory
		simout_dir2 (str): the other simOut directory
	Returns:
		diffs (dict): a dict describing the differences
	'''
	table1, table_message1 = open_table(simout_dir1, subdir)
	table2, table_message2 = open_table(simout_dir2, subdir)

	diffs = {}
	if table1 is None or table2 is None:
		diffs[subdir + '/'] = (table_message1, table_message2)
		return diffs

	# Compare Column names.
	column_names1 = set(table1.columnNames())
	column_names2 = set(table2.columnNames())

	# Diff the Table Columns.
	for key in column_names1 | column_names2:
		column1 = table1.readColumn2D(key) if key in column_names1 else Repr('<absent column>')
		column2 = table2.readColumn2D(key) if key in column_names2 else Repr('<absent column>')
		if is_repr(column1) or is_repr(column2):
			description = (column1, column2)
		else:
			description = compare_arrays(column1, column2)
		if description:
			diffs[subdir + '/' + key] = description

	# Diff the Table Attributes.
	attribute_names1 = set(table1.attributeNames())
	attribute_names2 = set(table2.attributeNames())

	for key in attribute_names1 | attribute_names2:
		v1 = table1.readAttribute(key) if key in attribute_names1 else Repr('<absent value>')
		v2 = table2.readAttribute(key) if key in attribute_names2 else Repr('<absent value>')
		if v1 != v2:
			# Call str(key) to avoid the u'...' unicode marker clutter
			diffs[subdir + '@' + str(key)] = (elide(v1), elide(v2))

	return diffs

def diff_simout(simout_dir1, simout_dir2):
	'''Diff two simOut dirs. Return a dict describing the differences.

	TODO(jerry): One could call this function in a Python shell and explore the
	values, but to do that it should return all the original values and make
	the caller elide them before printing.
	'''
	diffs = {}

	subdirs1 = set(os.listdir(simout_dir1))
	subdirs2 = set(os.listdir(simout_dir2))

	# ------------------------------------------------------------------------
	# Ignore the EvaluationTime table since CPU timing measurements will always
	# vary. It isn't really simulation output.
	#
	# Ignore the Daughter* cell inherited state pickle files. They don't
	# contain Tables and they only depend on wholecell/sim/divide_cell.py, not
	# the actual simulation.
	# ------------------------------------------------------------------------
	subdirs = (subdirs1 | subdirs2) - {'EvaluationTime'}

	for subdir in subdirs:
		if subdir.endswith('.cPickle'):
			continue

		diffs.update(diff_subdirs(subdir, simout_dir1, simout_dir2))

	return diffs

def cmd_diff_simout(simout_dir1, simout_dir2):
	'''Command line diff simout_dir2 against reference simout_dir1, with
	messages and returning an exit status code.
	'''
	if simout_dir1 == simout_dir2:
		return "diff_simouts: Don't diff a simOut directory against itself"

	print('Comparing: {}\n'.format((simout_dir1, simout_dir2)))

	if not op.isdir(simout_dir1):
		return 'diff_simouts: No simOut directory 1: {}'.format(simout_dir1)
	if not op.isdir(simout_dir2):
		return 'diff_simouts: No simOut directory 2: {}'.format(simout_dir2)

	diffs = diff_simout(simout_dir1, simout_dir2)

	if diffs:
		pprint(diffs)
		return 1
	else:
		print('The simOut dirs match')
		return 0


if __name__ == '__main__':
	if len(sys.argv) == 2 and sys.argv[1] in {'-h', '--help'}:
		print('''Usage:  diff_simouts <simOut1> <simOut2>
Diff two Whole Cell Model simOut directories

The output is a dict with keys like:

	* 'Main/': a subdir; present if there's a message about absent subdirs or
	  subdirs that don't contain Tables
	* 'Main@startTime': a Table attribute
	* 'RnaDegradationListener/DiffRelativeFirstOrderDecay': a Table column

The dict values are tuples of the corresponding simOut values or special
messages like:

	* <absent subdir>: an absent subdirectory
	* <absent column>: an absent column
	* <absent value>: an absent attribute value
	* Arrays are not equal (mismatch 21.6859279402%)...: A NumPy Testing
	  message summarizing the differences between two arrays
''')
		sys.exit()

	if len(sys.argv) < 3:
		sys.exit('diff_simouts: Need 2 simOut directory arguments')

	exit_status = cmd_diff_simout(sys.argv[1], sys.argv[2])
	sys.exit(exit_status)
