#! /usr/bin/env python
"""
Compare the output of two parca runs and show any differences.

Usage (SIMDIR is a path to simulation output containing a kb/ directory):
	runscripts/debug/compareParca.py SIMDIR1 SIMDIR2
"""

from __future__ import absolute_import, division, print_function

import argparse
import sys

from runscripts.reflect.object_tree import diff_trees, load_fit_tree, pprint_diffs


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Compare the Parca output between the two out/ sim-dirs and"
					" print a summary and count of the differences.")
	parser.add_argument('-c', '--count', action='store_true',
		help="Print just the diff line count, skipping the detailed diff lines.")
	parser.add_argument('sim_dir', nargs=2,
		help="The two out/ sim-dirs to compare.")

	args = parser.parse_args()
	dir1, dir2 = args.sim_dir

	tree1 = load_fit_tree(dir1)
	tree2 = load_fit_tree(dir2)

	diffs = diff_trees(tree1, tree2)
	pprint_diffs(diffs, print_diff_lines=not args.count)

	sys.exit(3 if diffs else 0)
