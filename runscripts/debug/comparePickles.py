#! /usr/bin/env python
"""
Compare two .cPickle files or all the .cPickle files in a pair of directories.
Show the differences or optionally just a count of difference lines.

Usage (PATH is a path like 'out/manual/intermediates'):
	runscripts/debug/comparePickles.py PATH1 PATH2
"""

import argparse
import os
import sys

from runscripts.reflect.object_tree import diff_dirs, diff_files
from wholecell.utils import constants


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Compare two .cPickle files"
					" or all the .cPickle files in two directories (in"
					" modification-time order)."
					" Print a count and optionally a summary of the differences.")
	parser.add_argument('-c', '--count', action='store_true',
		help="Print just the diff line count for each file, skipping the"
			 " detailed diff lines.")
	parser.add_argument('-f', '--final-sim-data', action='store_true',
		help="Append /kb/simData.cPickle to the two PATH args to make it a"
			 " little easier compare the final Parca output sim_data.")
	parser.add_argument('path', metavar='PATH', nargs=2,
		help="The two pickle files or directories to compare.")

	args = parser.parse_args()
	path1, path2 = args.path

	if args.final_sim_data:
		path1 = os.path.join(path1, constants.KB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)
		path2 = os.path.join(path2, constants.KB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)

	if os.path.isfile(path1):
		diff_count = diff_files(path1, path2, print_diff_lines=not args.count)
	else:
		diff_count = diff_dirs(path1, path2, print_diff_lines=not args.count)

	sys.exit(3 if diff_count else 0)
