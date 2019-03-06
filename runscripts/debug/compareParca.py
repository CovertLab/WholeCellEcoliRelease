from __future__ import absolute_import, division, print_function

import os
import cPickle
from pprint import pprint
import sys

from runscripts.reflect.object_tree import object_tree, diff_trees


def load_fit_tree(out_subdir):
	'''Load the Parca's output as an object_tree.'''
	path = os.path.join(
		os.getcwd(),
		'out',
		out_subdir,
		'kb',
		'simData_Fit_1.cPickle')

	with open(path, "rb") as f:
		sim_data = cPickle.load(f)

	return object_tree(sim_data)


# Compare the output of two parca runs, optionally running the parca first.

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('''Usage: {} SIMDIR1 SIMDIR2
Compare the Parca output between the two out/ sim-dirs and print a summary of
the differences.'''.format(sys.argv[0]))
		exit(2)

	dir1 = sys.argv[1]
	dir2 = sys.argv[2]

	once = load_fit_tree(dir1)
	twice = load_fit_tree(dir2)

	diffs = diff_trees(once, twice)
	pprint(diffs, width=160)

	exit(3 if diffs else 0)
