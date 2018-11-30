from __future__ import absolute_import, division, print_function

import os
import cPickle
from pprint import pprint

from runscripts.manual.runFitter import RunFitter
from runscripts.reflect.object_tree import object_tree, diff_trees

# Run the fitter and compare the new fitter output to the previously stored
# fitter output

if __name__ == '__main__':
	sim_data_file = os.path.join(
		os.getcwd(),
		'out',
		'manual',
		'kb',
		'simData_Fit_1.cPickle')

	with open(sim_data_file, "rb") as f:
		once = object_tree(cPickle.load(f))

	second = RunFitter()
	second.cli()

	with open(sim_data_file, "rb") as f:
		twice = object_tree(cPickle.load(f))

	recent = diff_trees(once, twice)
	newer = diff_trees(twice, once)

	pprint(recent)
	pprint(newer)
