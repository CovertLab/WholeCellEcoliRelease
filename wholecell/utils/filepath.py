"""
filepath.py
File and filename path utilities.
"""

from __future__ import absolute_import
from __future__ import division

import errno
import os

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
