"""
filepath.py
File and filename path utilities.
"""

from __future__ import absolute_import, division, print_function

import datetime
import errno
import json
import io
import os
import sys
from six.moves import range
if os.name == 'posix' and sys.version_info[0] < 3:
	# noinspection PyPackageRequirements
	import subprocess32 as subprocess2
	subprocess = subprocess2
else:
	import subprocess as subprocess3
	subprocess = subprocess3
from typing import Any, Generator, Optional, Sequence, Tuple

import six

import wholecell
from wholecell.utils.py3 import String


TIMEOUT = 60  # seconds

# The wcEcoli/ project root path which contains wholecell/.
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.realpath(wholecell.__file__)))
OUT_DIR = os.path.join(ROOT_PATH, 'out')
DEBUG_OUT_DIR = os.path.join(OUT_DIR, 'debug')

MATPLOTLIBRC_FILE = os.path.join(ROOT_PATH, 'matplotlibrc')

# Regex for current and previous timestamp() formats: 'YYYYMMDD.HHMMSS[.uuuuuu]'.
TIMESTAMP_PATTERN = r'\d{8}\.\d{6}(?:\.\d{6})?'

def makedirs(path, *paths):
	# type: (str, str) -> str
	"""Join one or more path components, make that directory path (using the
	default mode 0o0777), and return the full path.

	Raise FileExistsError if there's a file (not a directory) with that path.
	No exception if the directory already exists.
	"""
	full_path = os.path.join(path, *paths)

	if full_path:
		os.makedirs(full_path, exist_ok=True)

	return full_path

def timestamp(dt=None):
	# type: (Optional[datetime.datetime]) -> str
	"""Construct a datetime-timestamp from `dt` [default = `now()`], such as
	we use to timestamp a simulation output directory.
	"""
	if not dt:
		dt = datetime.datetime.now()

	return dt.strftime('%Y%m%d.%H%M%S')

def verify_file_exists(file_path, message=''):
	# type: (str, str) -> None
	"""Raise an IOError if file_path isn't an existing file."""
	if not os.path.isfile(file_path):
		raise IOError(errno.ENOENT,
			'Missing file "{}".  {}'.format(file_path, message))

def verify_dir_exists(dir_path, message=''):
	# type: (str, str) -> None
	"""Raise an IOError if dir_path isn't an existing directory."""
	if not os.path.isdir(dir_path):
		raise IOError(errno.ENOENT,
			'Missing dir "{}".  {}'.format(dir_path, message))

def run_cmd2(tokens, trim=True, timeout=TIMEOUT, env=None):
	# type: (Sequence[str], bool, Optional[int], Optional[dict]) -> Tuple[str, str]
	"""Run a shell command-line (in token list form) and return a tuple
	containing its (stdout, stderr).
	This does not expand filename patterns or environment variables or do other
	shell processing steps.

	Args:
		tokens: The command line as a list of string tokens.
		trim: Whether to trim off trailing whitespace. This is useful
			because the outputs usually end with a newline.
		timeout: timeout in seconds; None for no timeout.
		env: optional environment variables for the new process to use instead
			of inheriting the current process' environment.
	Returns:
		The command's stdout and stderr strings.
	Raises:
		OSError (e.g. FileNotFoundError [Python 3] or PermissionError),
		  subprocess.SubprocessError (TimeoutExpired or CalledProcessError)
	"""
	out = subprocess.run(
		tokens,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		check=True,
		env=env,
		encoding='utf-8',
		timeout=timeout)
	if trim:
		return out.stdout.rstrip(), out.stderr.rstrip()
	return out.stdout, out.stderr


def run_cmd(tokens, trim=True, timeout=TIMEOUT, env=None):
	# type: (Sequence[str], bool, Optional[int], Optional[dict]) -> str
	"""Run a shell command-line (in token list form) and return its stdout.
	See run_cmd2().
	"""
	return run_cmd2(tokens, trim=trim, timeout=timeout, env=env)[0]


def run_cmdline(line, trim=True, timeout=TIMEOUT, fallback=None):
	# type: (str, bool, Optional[int], Optional[str]) -> Optional[str]
	"""Run a shell command-line string then return its output or fallback if it
	failed. This does not expand filename patterns or environment variables or
	do other shell processing steps like quoting.

	Args:
		line: The command line as a string to split.
		trim: Whether to trim off trailing whitespace. This is useful
			because the subprocess output usually ends with a newline.
		timeout: timeout in seconds; None for no timeout.
		fallback: Return this if the subprocess fails, e.g. trying to run git
			in a Docker Image that has no git repo.
	Returns:
		The command's output string, or None if it couldn't even run.
	"""
	try:
		return run_cmd(tokens=line.split(), trim=trim, timeout=timeout)
	except (OSError, subprocess.SubprocessError) as e:
		if fallback is None:
			print('failed to run command line {}: {}'.format(line, e))
		return fallback


def git_hash():
	"""Return the source code git hash, or the environment variable
	$IMAGE_GIT_HASH if there's no git repo (in a Docker Image), or else '--'.
	"""
	return run_cmdline("git rev-parse HEAD",
					   fallback=os.environ.get("IMAGE_GIT_HASH", '--'))


def git_branch():
	"""Return the source code git branch name, or the environment variable
	$IMAGE_GIT_BRANCH if there's no git repo (in a Docker Image), or '--'.
	"""
	return run_cmdline("git symbolic-ref --short HEAD",
					   fallback=os.environ.get("IMAGE_GIT_BRANCH", '--'))

def write_file(filename, content):
	# type: (str, String) -> None
	"""Write text string `content` as a utf-8 text file."""
	with io.open(filename, 'w', encoding='utf-8') as f:
		f.write(six.text_type(content))

def write_json_file(filename, obj, indent=4):
	# type: (str, Any, int) -> None
	"""Write `obj` to a file in a pretty JSON format. This supports Unicode."""
	# Indentation puts a newline after each ',' so suppress the space there.
	message = json.dumps(obj, ensure_ascii=False, indent=indent,
		separators=(',', ': '), sort_keys=True) + '\n'
	write_file(filename, message)

def read_json_file(filename):
	# type: (str) -> Any
	"""Read and parse JSON file. This supports Unicode."""
	with io.open(filename, encoding='utf-8') as f:
		return json.load(f)

def iter_variants(variant_type, first_index, last_index):
	# type: (str, int, int) -> Generator[Tuple[int, str], None, None]
	"""Generate Variant subdirs (index, name) over [first .. last] inclusive."""
	# TODO(jerry): Return a list instead of generating items?
	for i in range(first_index, last_index + 1):
		yield i, os.path.join('{}_{:06d}'.format(variant_type, i))
