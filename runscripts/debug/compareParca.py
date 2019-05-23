from __future__ import absolute_import, division, print_function

import os
import cPickle
from pprint import pprint
import sys

from wholecell.utils import constants
from runscripts.reflect.object_tree import object_tree, diff_trees


def load_fit_tree(out_subdir):
	'''Load the parameter calculator's (Parca's) output as an object_tree.'''
	# For convenience, optionally add the prefix 'out/'.
	if not os.path.isabs(out_subdir) and not os.path.isdir(out_subdir):
		out_subdir = os.path.join('out', out_subdir)

	path = os.path.join(
		out_subdir,
		'kb',
		constants.SERIALIZED_SIM_DATA_FILENAME)

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
