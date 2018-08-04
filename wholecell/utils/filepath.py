"""
filepath.py
File and filename path utilities.
"""

from __future__ import absolute_import
from __future__ import division

import errno
import os
import subprocess


def makedirs(path, *paths):
	"""Join one or more path components, make that directory path (using the
	default mode 0o0777), and return the full path.

	Raise OSError if it can't achieve the result (e.g. the containing directory
	is readonly or the path contains a file); not if the directory already
	exists.
	"""
	full_path = os.path.join(path, *paths)

	try:
		os.makedirs(full_path)
	except OSError as e:
		if e.errno != errno.EEXIST or not os.path.isdir(full_path):
			raise

	return full_path

def run_cmd(tokens=None, line=''):
	"""Run a command-line program and return its output. Pass in `tokens` as a
	list of string tokens or else `line` as a string to split into tokens.

	This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
	Sherlock needs the latter to find libcrypto.so to run `git`.
	"""
	if not tokens:
		tokens = line.split()
	environ = {
		"PATH": os.environ["PATH"],
		"LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
		}
	out = subprocess.Popen(tokens, stdout = subprocess.PIPE, env=environ).communicate()[0]
	return out

def write_file(filename, content):
	"""Write string `content` as a text file."""
	with open(filename, "w") as f:
		f.write(content)
