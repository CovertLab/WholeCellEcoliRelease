from __future__ import absolute_import, division, print_function

import os
import cPickle
from pprint import pprint
import sys

from runscripts.manual.runFitter import RunFitter
from runscripts.reflect.object_tree import object_tree, diff_trees


def load_fit_tree(out_subdir):
	'''Load the Fitter's output as an object_tree.'''
	path = os.path.join(
		os.getcwd(),
		'out',
		out_subdir,
		'kb',
		'simData_Fit_1.cPickle')

	with open(path, "rb") as f:
		sim_data = cPickle.load(f)

	return object_tree(sim_data)


# Compare the output of two fitter runs, optionally running the fitter first.

if __name__ == '__main__':
	if len(sys.argv) > 2:    # diff the two named sim dirs
		dir1 = sys.argv[1]
		dir2 = sys.argv[2]
	elif len(sys.argv) > 1:  # diff the named sim dir against 'manual'
		dir1 = 'manual'
		dir2 = sys.argv[1]
	else:                    # run in 'new' then diff against 'manual'
		dir1 = 'manual'
		dir2 = 'new'

	once = load_fit_tree(dir1)

	if dir2 == 'new':  # TODO(jerry): Drop this RunFitter() feature?
		sys.argv[1:] = [dir2]
		second = RunFitter()
		second.cli()

	twice = load_fit_tree(dir2)

	recent = diff_trees(once, twice)
	pprint(recent)
