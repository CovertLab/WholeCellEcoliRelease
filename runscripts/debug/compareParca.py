#! /usr/bin/env python
"""
Compare the output of two parca runs and show any differences.

Usage (SIMDIR is a path to simulation output containing a kb/ directory):
	runscripts/debug/compareParca.py SIMDIR1 SIMDIR2
"""

from __future__ import absolute_import, division, print_function

from pprint import pprint
import sys

from runscripts.reflect.object_tree import diff_trees, load_fit_tree


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
