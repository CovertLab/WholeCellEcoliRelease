"""
Common code for scripts that manually run simulation and analysis operations
from the command line (e.g. outside of Fireworks workflows). This has code for
finding the simulation dir, variant subdirs, etc.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

import abc
import argparse
import datetime
import errno
import itertools
import re
import os
import pprint as pp
import time
import traceback
from typing import Any, Callable, Iterable, List, Optional, Tuple

import wholecell.utils.filepath as fp
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils.py3 import monotonic_seconds, process_time_seconds


METADATA_KEYS = (
	'timeline',
	'generations',
	'seed',
	'init_sims',
	'mass_distribution',
	'growth_rate_noise',
	'd_period_division',
	'variable_elongation_transcription',
	'variable_elongation_translation',
	'translation_supply',
	'trna_charging',
	'ppgpp_regulation',
	'superhelical_density',
	'recycle_stalled_elongation',
	'mechanistic_replisome',
	'mechanistic_aa_supply',
	'trna_attenuation',
	)

PARCA_KEYS = (
	'ribosome_fitting',
	'rnapoly_fitting',
	'cpus',
	'variable_elongation_transcription',
	'variable_elongation_translation')

SIM_KEYS = (
	'timeline',
	'length_sec',
	'timestep_safety_frac',
	'timestep_max',
	'timestep_update_freq',
	'log_to_disk_every',
	'jit',
	'mass_distribution',
	'growth_rate_noise',
	'd_period_division',
	'variable_elongation_transcription',
	'variable_elongation_translation',
	'translation_supply',
	'trna_charging',
	'ppgpp_regulation',
	'superhelical_density',
	'recycle_stalled_elongation',
	'mechanistic_replisome',
	'mechanistic_aa_supply',
	'trna_attenuation',
	'raise_on_time_limit',
	'log_to_shell',
	)

ANALYSIS_KEYS = (
	'plot',
	'cpus',
	'compile')

# Mapping of args that take a range and the args they will replace
RANGE_ARGS = {
	'variant_range': 'variant_index',
	'generation_range': 'generation',
	'seed_range': 'seed',
	# TODO: add daughter_range? - might require additional logic with increasing number per generation
	}

DEFAULT_VARIANT = ['wildtype', '0', '0']


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

def dashize(underscore):
	# type: (str) -> str
	return re.sub(r'_+', r'-', underscore)


class ScriptBase(metaclass=abc.ABCMeta):
	"""Abstract base class for scripts. This defines a template where
	`description()` describes the script,
	`define_parameters()` defines its command line parameters,
	`parse_args()` parses the command line args,
	`run()` does the work,
	`cli()` is the driving Command-Line Interpreter.
	"""

	# Regex to match a variant directory name. In the resulting match
	# object, group 1 is the variant_type and group 2 is the variant_index.
	VARIANT_DIR_PATTERN = re.compile(r'(.+)_(\d+)\Z')

	def __init__(self):
		self.range_options = []

	def description(self):
		# type: () -> str
		"""Describe the command line program. This defaults to the class name."""
		return type(self).__name__

	def epilog(self):
		# type: () -> Optional[str]
		"""Return the epilog help text; None for none."""
		return None

	def help(self):
		# type: () -> str
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
		flag if it isn't already defined. Overrides should call super.

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
		try:
			parser.add_argument('--verbose', action='store_true',
				help='Enable verbose logging.')
		except argparse.ArgumentError:
			pass  # ignore the conflict

	def define_parameter_bool(self, parser, underscore_name, default=False,
			help='', default_key=None):
		# type: (argparse.ArgumentParser, str, Any, str, Optional[str]) -> None
		"""Add a boolean option parameter to the parser. The CLI input can be
		`--name`, or its inverse `--no_name`. The default can be True or False, and
		changing it won't affect any of those explicit input forms. This method
		adds the default value to the help text.

		The parameter's `underscore_name` can contain underscores for easy
		searching in the code. This converts them to dashes for CLI convention.
		ArgumentParser does the reverse conversion when storing into `args`.

		Setting `default_key` copies the default value from
		`DEFAULT_SIMULATION_KWARGS[default_key]`.
		"""
		name = dashize(underscore_name)
		default = bool(
			DEFAULT_SIMULATION_KWARGS[default_key] if default_key else default)
		group = parser.add_mutually_exclusive_group()
		group.add_argument(
			'--' + name,
			default=default,
			action='store_true',
			help='Default = {}. {}'.format(default, help))
		group.add_argument(
			'--no-' + name,
			dest=underscore_name,
			action='store_false',
			help='Sets --{} to False'.format(name))

	def define_option(self, parser, underscore_name, datatype, default=None,
			help='', default_key=None, flag=''):
		# type: (argparse.ArgumentParser, str, Callable, Any, str, Optional[str], str) -> None
		"""Add an option with the given name and datatype to the parser.

		The parameter's `underscore_name` can contain underscores for easy
		searching in the code. This converts them to dashes for CLI convention.
		ArgumentParser does the reverse conversion when storing into `args`.

		Setting `default_key` copies the default value from
		`DEFAULT_SIMULATION_KWARGS[default_key]`.
		"""
		names = (('-' + flag,) if flag else ()) + (dashize('--' + underscore_name),)
		default = DEFAULT_SIMULATION_KWARGS[default_key] if default_key else default
		parser.add_argument(*names,
			type=datatype,
			default=default,
			help='({}; default {!r}) {}'.format(datatype.__name__, default, help)
			)

	def define_parameter_sim_dir(self, parser, default=None):
		# type: (argparse.ArgumentParser, Optional[Any]) -> None
		"""Add a `sim_dir` parameter to the command line parser. parse_args()
		will then use `args.sim_dir` to add `args.sim_path`.

		sim_dir identifies the simulation's output directory, defaulting to the
		latest timestamped (or else alphabetically first) subdirectory of
		"wcEcoli/out/".

		Call this in overridden define_parameters() methods as needed.
		"""
		parser.add_argument('sim_dir', nargs='?',
			default=default,
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
		parser.add_argument('-v', '--variant-index', type=int1,
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

	def define_elongation_options(self, parser):
		# type: (argparse.ArgumentParser) -> None
		"""Define the variable-elongation options for both Parca and Sim."""

		def add_bool_option(name, key, help):
			self.define_parameter_bool(parser, name, help=help, default_key=key)

		add_bool_option('variable_elongation_transcription', 'variable_elongation_transcription',
			help='Use a different elongation rate for different transcripts'
				 ' (currently increases rates for rRNA).'
				 ' Usually set this consistently between runParca and runSim.')
		add_bool_option('variable_elongation_translation', 'variable_elongation_translation',
			help='Use a different elongation rate for different polypeptides'
				 ' (currently increases rates for ribosomal proteins).'
				 ' Usually set this consistently between runParca and runSim.')

	def define_parca_options(self, parser, run_parca_option=False):
		# type: (argparse.ArgumentParser, bool) -> None
		"""Define Parca task options EXCEPT the elongation options."""

		if run_parca_option:
			self.define_parameter_bool(parser, 'run_parca', True,
				help='Run the Parca. The alternative, --no-run-parca, is useful'
					 ' to run more cell sims without rerunning the Parca.'
					 ' For that to work, the CLI args must specify the'
					 ' --timestamp and the same --description, --id, and'
					 ' --storage-root as a previous workflow that ran the Parca'
					 ' in order to locate its storage path. --no-run-parca makes'
					 ' other Parca CLI options irrelevant (the options below,'
					 ' through --no-debug-parca).')

		self.define_parameter_bool(parser, 'ribosome_fitting', True,
			help="Fit ribosome expression to protein synthesis demands.")
		self.define_parameter_bool(parser, 'rnapoly_fitting', True,
			help="Fit RNA polymerase expression to protein synthesis demands.")

		self.define_parameter_bool(parser, 'debug_parca', False,
			help='Make Parca calculate only one arbitrarily-chosen transcription'
				 ' factor condition when adjusting gene expression levels, leaving'
				 ' the other TFs at their input levels for faster Parca debugging.'
				 ' DO NOT USE THIS FOR A MEANINGFUL SIMULATION.')

	def define_sim_loop_options(self, parser, manual_script=False):
		# type: (argparse.ArgumentParser, bool) -> None
		"""Define options for running a series of sims."""

		# Variant
		parser.add_argument('-v', '--variant', nargs=3, default=DEFAULT_VARIANT,
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='''The variant type name, first index, and last index.
				See models/ecoli/sim/variants/__init__.py for the variant
				type choices and their supported index ranges, e.g.: wildtype,
				condition, meneParams, metabolism_kinetic_objective_weight,
				nutrientTimeSeries, and param_sensitivity.
				The meaning of the index values depends on the variant type. With
				wildtype, every index does the same thing, so it's a way to test
				that the simulation is repeatable.
				Default = ''' + ' '.join(DEFAULT_VARIANT))
		if manual_script:
			self.define_parameter_bool(parser, 'require_variants', False,
				help='''true => require the sim_data variant(s) specified by the
					--variant option to already exist; false => make the variant(s).
					Run makeVariants.py to make sim_data variants.''')

		# Simulation
		parser.add_argument('-g', '--generations', type=int, default=1,
			help='Number of cell sim generations to run per variant. (Single'
				 ' daughters only.) Default = 1')
		if manual_script:
			parser.add_argument(dashize('--total_gens'), type=int,
				help='(int) Total number of generations to write into the'
					 ' metadata.json file. Default = the value of --generations.')
		parser.add_argument('-s', '--seed', type=int, default=0,
			help="Simulation seed for the first generation of the first cell"
				 " lineage of every variant. The lineages (--init-sims) get"
				 " sequentially increasing seed numbers. The generations"
				 " (--generations) get seeds computed from the lineage seed and"
				 " the generation number. Default = 0")

		self.define_option(parser, 'init_sims', int, 1, flag='i',
			help='Number of initial sims (cell lineages) per variant. The'
				 ' lineages get sequential seeds starting with the --seed value.')

	def define_sim_options(self, parser):
		# type: (argparse.ArgumentParser) -> None
		"""Define sim task options EXCEPT the elongation options, with defaults
		from the sim class.
		"""

		def add_option(name, key, datatype, help):
			self.define_option(parser, name, datatype, help=help, default_key=key)

		def add_bool_option(name, key, help):
			self.define_parameter_bool(parser, name, help=help, default_key=key)

		self.define_option(parser, 'timeline', str, flag='t',
			default_key='timeline',
			help='The media timeline. See wholecell/utils/make_media.py,'
				 ' make_timeline() for timeline formatting details')
		add_option('length_sec', 'lengthSec', int,
			help='The maximum simulation time, in seconds. Useful for short'
				 ' simulations; not so useful for multiple generations.'
				 ' Default is 3 hours')
		add_option('timestep_safety_frac', 'timeStepSafetyFraction', float,
			help='Scale the time step by this factor if conditions are'
				 ' favorable, up the the limit of the max time step')
		add_option('timestep_max', 'maxTimeStep', float,
			help='the maximum time step, in seconds')
		add_option('timestep_update_freq', 'updateTimeStepFreq', int,
			help='frequency at which the time step is updated')
		add_option('log_to_disk_every', 'logToDiskEvery', int,
			help='frequency at which sim outputs are written to disk')

		add_bool_option('jit', 'jit',
			help='If true, jit compiled functions are used for certain'
				 ' processes, otherwise only uses lambda functions')
		add_bool_option('mass_distribution', 'massDistribution',
			help='If true, a mass coefficient is drawn from a normal distribution'
				 ' centered on 1; otherwise it is set equal to 1')
		add_bool_option('growth_rate_noise', 'growthRateNoise',
			help='If true, a growth rate coefficient is drawn from a normal'
				 ' distribution centered on 1; otherwise it is set equal to 1')
		add_bool_option('d_period_division', 'dPeriodDivision',
			help='If true, ends simulation once D period has occurred after'
				 ' chromosome termination; otherwise simulation terminates once'
				 ' a given mass has been added to the cell')
		add_bool_option('translation_supply', 'translationSupply',
			help='If true, the ribosome elongation rate is limited by the'
				 ' condition specific rate of amino acid supply; otherwise the'
				 ' elongation rate is set by condition')
		add_bool_option('trna_charging', 'trna_charging',
			help='if true, tRNA charging reactions are modeled and the ribosome'
				 ' elongation rate is set by the amount of charged tRNA	present.'
				 ' This option will override TRANSLATION_SUPPLY in the simulation.')
		add_bool_option('ppgpp_regulation', 'ppgpp_regulation',
			help='if true, ppGpp concentration is determined with kinetic equations.')
		add_bool_option('superhelical_density', 'superhelical_density',
			help='if true, dynamically calculate superhelical densities of each DNA segment')
		add_bool_option('recycle_stalled_elongation', 'recycle_stalled_elongation',
						help='if true, recycle RNAP and fragment bases when transcription'
							 'elongation is stalled in ntp-limiting conditions')
		add_bool_option('mechanistic_replisome', 'mechanistic_replisome',
			help='if true, replisome initiation is mechanistic (requires'
				 ' appropriate number of subunits to initiate)')
		add_bool_option('mechanistic_aa_supply', 'mechanistic_aa_supply',
			help='if true, amino acid supply is mechanistic (depends on'
				 ' concentrations of enzymes and amino acids)')
		add_bool_option('trna_attenuation', 'trna_attenuation',
			help='if true, transcriptional attenuation by charged tRNA is enabled')
		add_bool_option('raise_on_time_limit', 'raise_on_time_limit',
			help='if true, the simulation raises an error if the time limit'
				 ' (--length-sec) is reached before division.')
		add_bool_option('log_to_shell', 'logToShell',
			help='if true, logs output to the shell')

	def define_range_options(self, parser, *range_keys):
		# type: (argparse.ArgumentParser, *str) -> None
		"""Adds options to the arg parser for values that can take on ranges."""

		for key in range_keys:
			option = '{}_range'.format(key)
			upper = key.upper()
			override = dashize('--{}'.format(RANGE_ARGS[option]))
			parser.add_argument(dashize('--{}'.format(option)), nargs=2, default=None, type=int,
				metavar=('START_{}'.format(upper), 'END_{}'.format(upper)),
				help='The range of variants to run.  Will override {} option.'.format(override))
			self.range_options.append(option)

	def parse_args(self):
		# type: () -> argparse.Namespace
		"""Parse the command line args: Construct an ArgumentParser, call
		`define_parameters()` to define parameters including subclass-specific
		parameters, use it to parse the command line into an
		`argparse.Namespace`, and return that.

		(A `Namespace` is an object with attributes and some methods like
		`__repr__()` and `__eq__()`. Call `vars(args)` to turn it into a dict.)
		"""
		parser = argparse.ArgumentParser(
			description=self.help(),
			epilog=self.epilog())

		self.define_parameters(parser)

		return parser.parse_args()

	def update_args(self, args):
		# type: (argparse.Namespace) -> None
		"""
		Update or add to parsed arguments.

		If there's a `sim_dir` arg [see define_parameter_sim_dir()], set
		`args.sim_path`. If there's also a `variant_index` arg [see
		define_parameter_variant_index()], set `args.variant_dir` to the
		find_variant_dir() tuple.

		When overriding, first call super().
		"""

		if 'sim_dir' in args:
			args.sim_path = find_sim_path(args.sim_dir)

			if 'variant_index' in args:
				args.variant_dir = self.find_variant_dir(
					args.sim_path, args.variant_index)

	def extract_range_args(self, args):
		# type: (argparse.Namespace) -> List[List[int]]
		"""
		Extracts arguments that have been specified as ranges for other arguments.

		Returns:
			lists of possible values that each argument with a range can take,
			ordered according to self.range_options
		"""

		range_args = []
		for range_option in self.range_options:
			if getattr(args, range_option):
				start, end = getattr(args, range_option)
				values = list(range(start, end+1))
			else:
				values = [getattr(args, RANGE_ARGS[range_option])]

			range_args.append(values)

		return range_args

	def set_range_args(self, args, params):
		# type: (argparse.Namespace, Iterable[int]) -> None
		"""Sets arguments from a combination of values from ranges."""

		for range_id, param in zip(self.range_options, params):
			setattr(args, RANGE_ARGS[range_id], param)

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

		exceptions = []
		range_args = self.extract_range_args(self.parse_args())
		# TODO (Travis): have option to parallelize when ranges given
		# TODO (Travis): check if variant/seed/gen range combo does not exist
		for params in itertools.product(*range_args):
			# Start with original args for each iteration since update_args
			# overwrites some args and might handle undefined values differently
			# with different values from ranges. There is not a good way to copy
			# the unmodified Namespace without using copy.deepcopy().
			args = self.parse_args()
			self.set_range_args(args, params)
			self.update_args(args)

			location = getattr(args, 'sim_path', '')
			if location:
				location = ' at ' + location

			start_real_sec = monotonic_seconds()
			print('{}: {}{}'.format(time.ctime(), self.description(), location))
			pp.pprint({'Arguments': vars(args)})

			start_process_sec = process_time_seconds()
			try:
				self.run(args)
			except Exception as e:
				# Handle exceptions after completion if running multiple params
				if range_args and max([len(a) for a in range_args]) > 1:
					traceback.print_exc()
					exceptions.append((params, e))
				else:
					raise
			elapsed_process = process_time_seconds() - start_process_sec

			elapsed_real_sec = monotonic_seconds() - start_real_sec
			print("{}: Elapsed time {:1.2f} sec ({}); CPU {:1.2f} sec".format(
				time.ctime(),
				elapsed_real_sec,
				datetime.timedelta(seconds=elapsed_real_sec),
				elapsed_process,
				))

		# Handle any exceptions that occurred
		if exceptions:
			for params, ex in exceptions:
				param_str = ', '.join(['{}: {}'.format(RANGE_ARGS[o], p)
					for o, p in zip(self.range_options, params)])
				print('Error with param set ({}): "{}"'.format(param_str, ex))

			raise RuntimeError('Exception in one or more parameter sets (see above).')


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
