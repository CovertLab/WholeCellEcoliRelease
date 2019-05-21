"""
Common code for scripts that manually run simulation and analysis operations
from the command line (e.g. outside of Fireworks workflows). This has code for
finding the simulation dir, variant subdirs, etc.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import abc
import argparse
import datetime
import errno
import re
import os
import pprint as pp
import time
from typing import Any, Callable, List, Optional, Tuple

import wholecell.utils.filepath as fp


def default_wcecoli_out_subdir_path():
	# type: () -> str
	"""Return an absolute path to the most interesting subdirectory of
	wcEcoli/out/: the subdirectory name that starts with the latest timestamp
	or (if none) the one that's first alphabetically.
	"""
	out_dir = os.path.join(fp.ROOT_PATH, 'out')
	timestamped = re.compile(fp.TIMESTAMP_PATTERN)
	fallback = None

	for directory in sorted(os.listdir(out_dir), reverse=True):
		path = os.path.join(out_dir, directory)
		if os.path.isdir(path):
			if timestamped.match(directory):
				return path
			fallback = path

	if fallback:
		return fallback

	raise IOError(errno.ENOENT,
		'"{}" has no subdirectories.  Run runParca?'.format(out_dir))

def find_sim_path(directory=None):
	# type: (Optional[str]) -> str
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
		input_dir = os.path.join(fp.ROOT_PATH, directory)
	else:
		input_dir = os.path.join(fp.ROOT_PATH, 'out', directory)

	fp.verify_dir_exists(input_dir, "Need a simulation dir.")
	return input_dir

def str_to_bool(s):
	# type: (str) -> bool
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
	VARIANT_DIR_PATTERN = re.compile(r'(.+)_(\d+)\Z')

	def description(self):
		"""Describe the command line program. This defaults to the class name."""
		return type(self).__name__

	def help(self):
		"""Return help text for the Command Line Interface. This defaults to a
		string constructed around `self.description()`.
		"""
		return 'Run {}.'.format(self.description())

	def list_variant_dirs(self, sim_path):
		# type: (str) -> List[Tuple[str, str, int]]
		"""List the available variant subdirectories of the given sim_path,
		in alphabetical order,
		returning for each a tuple (subdir_name, variant_type, variant_index),
		where the variant_index is an int.
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
		# type: (argparse.ArgumentParser) -> None
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
		# type: (argparse.ArgumentParser, str, Any, str) -> None
		"""Add a boolean option parameter to the parser. The CLI input can be
		`--name`, `--no_name`, `--name true`, `--name false`, `--name 1`,
		`--name 0`, `--name=true`, etc. The default can be True or False, and
		changing it won't affect any of those explicit input forms. This method
		adds the default value to the help text.
		"""
		default = bool(default)
		examples = 'true or 1' if default else 'false or 0'
		group = parser.add_mutually_exclusive_group()
		group.add_argument('--' + name, nargs='?', default=default,
			const='true',  # needed for nargs='?'
			type=str_to_bool,
			help='({}; {}) {}'.format('bool', examples, help))
		group.add_argument('--no_' + name, dest=name, action='store_false',
			help='Like {}=0'.format(name))

	def define_option(self, parser, name, datatype, default, help):
		# type: (argparse.ArgumentParser, str, Callable, Any, str) -> None
		"""Add an option with the given name and datatype to the parser."""
		parser.add_argument('--' + name,
			type=datatype,
			default=default,
			help='({}; {}) {}'.format(datatype.__name__, default, help)
			)

	def define_parameter_sim_dir(self, parser):
		# type: (argparse.ArgumentParser) -> None
		"""Add a `sim_dir` parameter to the command line parser. parse_args()
		will then use `args.sim_dir` to add `args.sim_path`.

		sim_dir identifies the simulation's output directory, defaulting to the
		latest timestamped (or else alphabetically first) subdirectory of
		"wcEcoli/out/".

		Call this in overridden define_parameters() methods as needed.
		"""
		parser.add_argument('sim_dir', nargs='?',
			help='''The simulation "out/" subdirectory to read from (optionally
				starting with "out/"), or an absolute directory name, or
				default to the "out/" subdirectory name that starts with
				the latest timestamp or else the one that's first
				alphabetically.''')

	def define_parameter_variant_index(self, parser):
		# type: (argparse.ArgumentParser) -> None
		"""Add a `variant_index` parameter to the command line parser.
		parse_args() will then use the `variant_index` and `sim_path`
		arguments, call find_variant_dir(), and set `args.variant_dir` to the
		first matching variant.

		Call this in overridden define_parameters() methods as needed.
		"""
		int1 = int  # type: Callable
		parser.add_argument('-v', '--variant_index', type=int1,
			help='The simulation variant number (int), e.g. 1 to find a'
				 ' subdirectory like "condition_000001".')

	def find_variant_dir(self, sim_path, index=None):
		# type: (str, Optional[int]) -> Tuple[str, str, int]
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
		# type: () -> argparse.Namespace
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
		# type: (argparse.Namespace) -> None
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

		start_wall_sec = time.time()
		print('{}: {}{}'.format(
			time.ctime(start_wall_sec), self.description(), location))
		pp.pprint({'Arguments': vars(args)})

		start_process_sec = time.clock()
		self.run(args)
		end_process_sec = time.clock()
		elapsed_process = end_process_sec - start_process_sec

		end_wall_sec = time.time()
		elapsed_wall = end_wall_sec - start_wall_sec
		print("{}: Elapsed time {:1.2f} sec ({}); {:1.2f} sec in process".format(
			time.ctime(end_wall_sec),
			elapsed_wall,
			datetime.timedelta(seconds=elapsed_wall),
			elapsed_process,
			))


class TestScript(ScriptBase):
	"""To test out the command line parser."""

	def define_parameters(self, parser):
		super(TestScript, self).define_parameters(parser)
		parser.add_argument('--seed', default='000001', help='simulation seed')

	def run(self, args):
		print("[TEST] Run args:", args)


if __name__ == '__main__':
	script = TestScript()
	script.cli()
