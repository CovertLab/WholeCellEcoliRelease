#!/usr/bin/env python

"""
Summarize info on the Python runtime environment that's useful for
debugging and performance tuning.
"""

from __future__ import absolute_import, division, print_function

import multiprocessing
import numpy as np
import os
import scipy
import sys

from wholecell.utils import parallelization


def subtitle(module):
	"""Print a section subtitle with the module's name and version."""
	print()
	print(f"Module: {getattr(module, '__name__', '(no name)')}"
		  f" {getattr(module, '__version__', '')}")

def print_dictionary(name, a_dict, select_keys):
	"""Print select keys from a dict, as available, in order, one per line."""
	print("{}:".format(name))
	for key in select_keys:
		value = "'{}'".format(a_dict[key]) if key in a_dict else '--'
		print("  '{}': {}".format(key, value))

def main():
	subtitle(multiprocessing)
	print("multiprocessing.cpu_count(): {}".format(multiprocessing.cpu_count()))
	print(f'parallelization.cpus(): {parallelization.cpus()}')

	subtitle(np)
	np.show_config()

	subtitle(os)
	print_dictionary("os.environ", os.environ, [
		'OPENBLAS_NUM_THREADS', 'HOME', 'LIBRARY_PATH', 'PI_HOME', 'PYENV_ROOT',
		'PYTHONPATH', 'SHERLOCK'])
	print("os.getcwd(): {}".format(os.getcwd()))
	print("os.uname(): {}".format(os.uname()))

	subtitle(scipy)
	scipy.__config__.show()

	subtitle(sys)
	print("sys.platform: {}".format(sys.platform))
	print("sys.prefix: {}".format(sys.prefix))
	print("sys.version: {}".format(sys.version))
	print("sys.api_version: {}".format(sys.api_version))
	print("sys.version_info: {}".format(sys.version_info))
	print()


if __name__ == "__main__":
	main()
