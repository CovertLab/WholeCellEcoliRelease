#!/usr/bin/env python

"""
Summarize info on the Python runtime environment that's useful for
debugging and performance tuning.
"""

import multiprocessing
import os
import scipy
import sys

def main():
	print
	print "multiprocessing.cpu_count(): {}".format(multiprocessing.cpu_count())
	print

	select_env = {key: value for (key, value) in os.environ.iteritems()
				  if key in {'HOME', 'PYENV_ROOT'}}
	print "os.environ: {}...".format(select_env)
	print "os.getcwd(): {}".format(os.getcwd())
	print "os.uname(): {}".format(os.uname())
	print

	print "scipy.__config__.show()"
	scipy.__config__.show()
	print

	print "sys.platform(): {}".format(sys.platform)
	print "sys.prefix: {}".format(sys.prefix)
	print "sys.version: {}".format(sys.version)
	print "sys.api_version: {}".format(sys.api_version)
	print "sys.version_info: {}".format(sys.version_info)
	print


if __name__ == "__main__":
	main()
