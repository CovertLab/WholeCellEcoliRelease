"""
Common code for scripts that manually run simulation and analysis operations
from the command line (e.g. outside of Fireworks workflows). This has code for
finding the simulation dir, variant subdirs, etc.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import abc
import argparse
import datetime
import errno
import re
import os
import time

import wholecell


# The wcEcoli project root path.
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(wholecell.__file__)))


def default_wcecoli_out_subdir_path():
	"""Return an absolute path to the most interesting subdirectory of
	wcEcoli/out/: the one that starts with the latest timestamp or else the
	alphabetically first subdirectory.
	"""
	out_dir = os.path.join(ROOT_PATH, 'out')
	fallback = None

	for directory in sorted(os.listdir(out_dir), reverse=True):
		path = os.path.join(out_dir, directory)
		if os.path.isdir(path):
			if directory[0].isdigit():
				return path
			fallback = path

	if fallback:
		return fallback

	raise IOError(errno.ENOENT,
		'"{}" has no subdirectories.  Run the Fitter?'.format(out_dir))

def find_sim_path(directory=None):
	"""Find a simulation path, looking for the given directory name as an
	absolute path, or as a subdirectory of wcEcoli/out/, or as a subdirectory
	name that starts with out/, or (if None) call
	default_wcecoli_out_subdir_path().
	"""
	if directory is None:
		input_dir = default_wcecoli_out_subdir_path()
	elif os.path.isabs(directory):
		input_dir = directory
	elif directory.startswith('out/'):
		input_dir = os.path.join(ROOT_PATH, directory)
	else:
		input_dir = os.path.join(ROOT_PATH, 'out', directory)

	if not os.path.isdir(input_dir):
		raise IOError(errno.ENOENT, '{} is not a simulation path'.format(input_dir))
	return input_dir

def str_to_bool(s):
	"""Convert a string command line parameter value to a bool. This ignores
	case, accepting true, false, 1, or 0.
	"""
	s = s.lower()
	if s not in {'true', 'false', '1', '0'}:
		raise ValueError('Expected a bool, not %s' % s)
	return s in {'true', '1'}


class ScriptBase(object):
	"""Abstract base class for scripts. This defines a template where
	`description()` describes the script,
	`define_parameters()` defines its command line parameters,
	`parse_args()` parses the command line args,
	`run()` does the work,
	`cli()` is the driving Command-Line Interpreter.
	"""
	__metaclass__ = abc.ABCMeta

	# Regex to match a variant directory name. In the resulting match
	# object, group 1 is the variant_type and group 2 is the variant_index.
	VARIANT_DIR_PATTERN = re.compile(r'([a-zA-Z_\d]+)_(\d+)\Z')

	def description(self):
		"""Describe the command line program. This defaults to the class name."""
		return type(self).__name__

	def help(self):
		"""Return help text for the Command Line Interface. This defaults to a
		string constructed around `self.description()`.
		"""
		return 'Run {}.'.format(self.description())

	def timestamp(self, dt=None):
		"""Construct a datetime-timestamp from `dt`; default = now()."""
		if not dt:
			dt = datetime.datetime.now()

		# TODO: Simplify to `format(datetime_value, '%Y%m%d.%H%M%S.%f')`?
		return "%04d%02d%02d.%02d%02d%02d.%06d" % (
			dt.year, dt.month, dt.day,
			dt.hour, dt.minute, dt.second,
			dt.microsecond)

	def list_variant_dirs(self, sim_path):
		"""List the available variant subdirectories of the given sim_path,
		returning for each a tuple (subdir_name, variant_type, variant_index),
		with the variant_index as an int.
		"""
		available = []

		for subdir in sorted(os.listdir(sim_path)):
			match = self.VARIANT_DIR_PATTERN.match(subdir)
			if match:
				path = os.path.join(sim_path, subdir)
				if os.path.isdir(path):
					available.append((subdir, match.group(1), int(match.group(2))))

		return available

	def define_parameters(self, parser):
		"""Define command line parameters. This base method defines a --verbose
		flag. Overrides should call super.

		Examples include positional arguments
			`parser.add_argument('variant', nargs='?',
			help='Simulation variant.')`
		options
			`parser.add_argument('--seed', default='000000',
			help='Simulation seed.')`.
		and flags
			`parser.add_argument('--verbose', action='store_true',
			help='Enable verbose logging.')`.
		"""
		parser.add_argument('--verbose', action='store_true',
			help='Enable verbose logging.')

	def define_parameter_bool(self, parser, name, default, help):
		"""Add a boolean option parameter to the parser. The CLI input can be
		`--name`, `--no_name`, `--name true`, `--name false`, `--name 1`,
		`--name 0`, `--name=true`, etc. The default can be True or False, and
		changing it won't affect any of those explicit input forms. This method
		adds the default value to the help text.
		"""
		default = bool(default)
		group = parser.add_mutually_exclusive_group()
		group.add_argument('--' + name, nargs='?', default=default,
			const='true',  # needed for nargs='?'
			type=str_to_bool,
			help='{}. Default = {}'.format(help, default))
		group.add_argument('--no_' + name, dest=name, action='store_false')

	def define_parameter_sim_dir(self, parser):
		"""Add a `sim_dir` parameter to the command line parser. parse_args()
		will then use `args.sim_dir` to add `args.sim_path`.

		sim_dir identifies the simulation's output directory, defaulting to the
		latest timestamped (or else alphabetically first) subdirectory of
		"wcEcoli/out/".

		Call this in overridden define_parameters() methods as needed.
		"""
		parser.add_argument('sim_dir', nargs='?',
			help='The simulation "out/" subdirectory to read from (optionally'
				 ' starting with "out/"), or an absolute directory name, or'
				 ' default to the the most interesting subdirectory of "out/".')

	def define_parameter_variant_index(self, parser):
		"""Add a `variant_index` parameter to the command line parser.
		parse_args() will then use the `variant_index` and `sim_path`
		arguments, call find_variant_dir(), and set `args.variant_dir`.

		Call this in overridden define_parameters() methods as needed.
		"""
		parser.add_argument('-v', '--variant_index', type=int,
			help='The simulation variant number (int), e.g. 1 to find a'
				 ' subdirectory like "condition_000001".')

	def find_variant_dir(self, sim_path, index=None):
		"""Find a simulation variant dir in the given `sim_path` for the given
		`index`, returning a tuple (subdir_name, variant_type, variant_index)
		or raising IOError. If `index` is None, return the first available
		variant, otherwise return the first variant with the given index (int).
		"""
		available = self.list_variant_dirs(sim_path)

		if index is None and available:
			return available[0]

		for choice in available:
			if index == choice[2]:
				return choice

		raise IOError(errno.ENOENT,
			'Simulation variant directory not found for variant index {} in'
			' sim_path {}'.format(
				'(any)' if index is None else index,
				sim_path))

	def parse_args(self):
		"""Parse the command line args: Construct an ArgumentParser, call
		`define_parameters()` to define parameters including subclass-specific
		parameters, use it to parse the command line into an
		`argparse.Namespace`, and return that.

		If there's a `sim_dir` arg [see define_parameter_sim_dir()], set
		`args.sim_path`. If there's also a `variant_index` arg [see
		define_parameter_variant_index()], set `args.variant_dir` to the
		find_variant_dir() tuple.

		When overriding, first call super().

		(A `Namespace` is an object with attributes and some methods like
		`__repr__()` and `__eq__()`. Call `vars(args)` to turn it into a dict.)
		"""
		parser = argparse.ArgumentParser(description=self.help())

		self.define_parameters(parser)

		args = parser.parse_args()

		if 'sim_dir' in args:
			args.sim_path = find_sim_path(args.sim_dir)

			if 'variant_index' in args:
				args.variant_dir = self.find_variant_dir(
					args.sim_path, args.variant_index)

		return args

	@abc.abstractmethod
	def run(self, args):
		"""Run the operation with the given arguments. If args.verbose,
		overrides can do verbose logging.
		"""
		raise NotImplementedError("ScriptBase subclass must implement run()")

	def cli(self):
		"""Command Line Interpreter: parse_args() then run(). This also prints
		a starting message (including args.sim_path if defined) and an ending
		message (including the elapsed run time).
		"""
		args = self.parse_args()

		location = getattr(args, 'sim_path', '')
		if location:
			location = ' at ' + location

		print '{}: {}{}'.format(time.ctime(), self.description(), location)
		if args.verbose:
			print '    args: {}'.format(args)

		start_sec = time.clock()
		self.run(args)
		end_sec = time.clock()
		elapsed = datetime.timedelta(seconds = (end_sec - start_sec))

		print "Run in {}h {}m {}s total".format(*str(elapsed).split(':'))


class TestScript(ScriptBase):
	"""To test out the command line parser."""

	def define_parameters(self, parser):
		super(TestScript, self).define_parameters(parser)
		parser.add_argument('--seed', default='000001', help='simulation seed')

	def run(self, args):
		print "[TEST] Run args:", args


if __name__ == '__main__':
	script = TestScript()
	script.cli()
