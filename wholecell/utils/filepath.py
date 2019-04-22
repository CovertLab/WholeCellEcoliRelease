"""
filepath.py
File and filename path utilities.
"""

from __future__ import absolute_import
from __future__ import division

import datetime
import errno
import json
import io
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

def timestamp(dt=None):
	"""Construct a datetime-timestamp from `dt` [default = `now()`], such as
	we use to timestamp a simulation output directory.
	"""
	if not dt:
		dt = datetime.datetime.now()

	# TODO: Simplify to `format(datetime_value, '%Y%m%d.%H%M%S.%f')`?
	return "%04d%02d%02d.%02d%02d%02d.%06d" % (
		dt.year, dt.month, dt.day,
		dt.hour, dt.minute, dt.second,
		dt.microsecond)

def run_cmd(tokens, trim=True):
	"""Run a shell command-line (in token list form) and return its output.
	This does not expand filename patterns or environment variables or do other
	shell processing steps.

	This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
	Sherlock needs the latter to find libcrypto.so to run `git`.

	Args:
		tokens (list): The command line as a list of string tokens.
		trim (bool): Whether to trim off trailing whitespace. This is useful
			because the subprocess output usually ends with a newline.
	"""
	environ = {
		"PATH": os.environ["PATH"],
		"LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
		}
	out = subprocess.Popen(tokens, stdout = subprocess.PIPE, env=environ).communicate()[0]
	if trim:
		out = out.rstrip()
	return out

def run_cmdline(line, trim=True):
	"""Run a shell command-line string and return its output. This does not
	expand filename patterns or environment variables or do other shell
	processing steps.

	This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
	Sherlock needs the latter to find libcrypto.so to run `git`.

	Args:
		line (str): The command line as a string.
		trim (bool): Whether to trim off trailing whitespace. This is useful
			because the subprocess output usually ends with a newline.
	"""
	try:
		return run_cmd(tokens=line.split(), trim=trim)
	except StandardError as e:
		print('failed to run command line {}: {}'.format(line, e))
		return None

def write_file(filename, content):
	"""Write string `content` as a text file."""
	with open(filename, "w") as f:
		f.write(content)

def write_json_file(filename, obj, indent=4):
	"""Write `obj` to a file in a pretty JSON format. This supports Unicode."""
	# Indentation puts a newline after each ',' so suppress the space there.
	message = json.dumps(obj, ensure_ascii=False, indent=indent,
		separators=(',', ': '), sort_keys=True) + '\n'
	write_file(filename, message.encode('utf-8'))

def read_json_file(filename):
	"""Read and parse JSON file. This supports Unicode."""
	with io.open(filename, encoding='utf-8') as f:
		return json.load(f)
