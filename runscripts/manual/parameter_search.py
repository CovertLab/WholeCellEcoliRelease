#! /usr/bin/env python
"""
Run a parameter search using parca and simulation output.  See methods and
solvers in models/ecoli/sim/parameter_search and wholecell/optimization,
respectively.

TODO: Share more code with fw_queue.py and other runscripts.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from models.ecoli.sim.parameter_search import PARAMETER_METHODS
from wholecell.optimization import SOLVERS
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp


# Command line arg defaults for solver options
DEFAULT_ITERATIONS = 5
DEFAULT_STARTING_ITERATION = 0
DEFAULT_LEARNING_RATE = 1
DEFAULT_PARAMETER_STEP = 0.1
DEFAULT_MAX_CHANGE = 0.1
DEFAULT_ALPHA = 0.1
DEFAULT_GAMMA = 0.1


class RunParameterSearch(scriptBase.ScriptBase):
	"""Runs various parameter search algorithms to optimize desired objectives."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli parameter search'

	def help(self):
		"""Return help text for the Command Line Interface."""
		return '''Run a {}.
				If the sim_path ends with a dir like
				"20190704.101500__Latest_sim_run", this will get the
				timestamp and description from the path to write into
				metadata.json.
				The command line option names are long but you can use any
				unambiguous prefix.'''.format(self.description())

	def define_parameters(self, parser):
		super().define_parameters(parser)

		# TODO: reuse params from other parsers?
		# self.define_parameter_sim_dir(parser)
		# self.define_sim_loop_options(parser, manual_script=True)
		# self.define_sim_options(parser)
		# self.define_elongation_options(parser)

		default_solver = list(SOLVERS.keys())[0]

		self.define_parameter_sim_dir(parser)
		parser.add_argument('--timestamp', action='store_true',
			help='Timestamp the given `sim_dir`, transforming e.g.'
				 ' "Fast run" to "20190514.135600__Fast_run".')

		parser.add_argument('--solver',
			default=default_solver,
			choices=SOLVERS.keys(),
			help=f'Solver for optimizing parameters (default: {default_solver}).')
		parser.add_argument('--method',
			required=True,
			choices=PARAMETER_METHODS.keys(),
			help='Class defining parameters and an objective for parameter search.')
		self.define_option(parser, 'iterations', int, default=DEFAULT_ITERATIONS,
			help=f'Number of iterations to update parameters before stopping.')
		self.define_option(parser, 'starting_iteration', int, default=DEFAULT_STARTING_ITERATION,
			help=f'Iteration number to begin with (not guaranteed to reproduce results exactly).')
		self.define_option(parser, 'learning_rate', float, default=DEFAULT_LEARNING_RATE,
			help=f'Learning rate for updating parameters.')
		self.define_option(parser, 'parameter_step', float, default=DEFAULT_PARAMETER_STEP,
			help=f'Fraction to update parameters by to determine the gradient.')
		self.define_option(parser, 'max_change', float, default=DEFAULT_MAX_CHANGE,
			help=f'Maximum fraction to update a parameter in a given time step.')
		self.define_option(parser, 'cpus', int, default=1, flag='c',
			help=f'Number of CPUs to use for running parca/sims in parallel for a given iteration.')

		# SPSA options
		self.define_option(parser, 'alpha', float, default=DEFAULT_ALPHA,
			help=f'Set alpha parameter for SPSA solver to control learning rate decay.')
		self.define_option(parser, 'gamma', float, default=DEFAULT_GAMMA,
			help=f'Set gamma parameter for SPSA solver to control parameter update decay.')

	def parse_args(self):
		args = super().parse_args()

		if args.timestamp:
			args.sim_dir = fp.timestamp() + '__' + args.sim_dir.replace(' ', '_')

		args.sim_path = scriptBase.find_sim_path(directory=args.sim_dir, makedirs=True)

		return args

	def run(self, args):
		fp.makedirs(args.sim_path, constants.KB_DIR)
		fp.makedirs(args.sim_path, constants.METADATA_DIR)  # TODO: write metadata?

		method = PARAMETER_METHODS[args.method]()
		solver = SOLVERS[args.solver](method, args)

		# TODO: address this issue by saving raw/sim data directly or saving the updates to apply
		if args.starting_iteration != 0:
			print('Warning: starting at a new iteration is not guaranteed to give identical results.'
				' The updated parameters should be saved after each iteration to restart at the same'
				' place in the algorithm.')

		for it in range(args.starting_iteration, args.starting_iteration + args.iterations):
			print(f'** Starting iteration {it} **')
			objective = solver.run_iteration()
			solver.print_update(objective)


if __name__ == '__main__':
	script = RunParameterSearch()
	script.cli()
