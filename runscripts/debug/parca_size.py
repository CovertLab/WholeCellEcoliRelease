#! /usr/bin/env python
"""
Find the size of attributes in a sim_data object.  Shows totals of all
attributes contained within each attribute.

Usage (SIMDIR is a path to simulation output containing a kb/ directory):
	runscripts/debug/parca_size.py SIMDIR
"""

from __future__ import absolute_import, division, print_function

from pprint import pprint
import sys

from runscripts.reflect.object_tree import load_fit_tree, size_tree


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('''Usage: {} SIMDIR
Find size of attributes in a sim_data object'''.format(sys.argv[0]))
		exit(2)

	sim_data = load_fit_tree(sys.argv[1])
	sizes = size_tree(sim_data, 0.1)
	pprint(sizes, width=160)
