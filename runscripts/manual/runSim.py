"""
Run wcEcoli cell simulations, supporting multiple variants, multiple initial
seeds, and multiple generations but only single daughters per generation.

Prerequisite: Run the parameter calculator (runFitter.py).

TODO: Share lots of code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import cPickle
import re
import os
from typing import Tuple

from wholecell.fireworks.firetasks import SimulationDaughterTask, SimulationTask, VariantSimDataTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp


SIM_DIR_PATTERN = r'({})__(.+)'.format(fp.TIMESTAMP_PATTERN)


def parse_timestamp_description(sim_path):
	# type: (str) -> Tuple[str, str]
	"""Parse `timestamp, description` from a sim_path that ends with a dir like
	'20190704.101500__Latest_sim_run' or failing that, return defaults.
	"""
	sim_dir = os.path.basename(sim_path)
	if not sim_dir:  # sim_path is empty or ends with '/'
		sim_dir = os.path.basename(os.path.dirname(sim_path))

	match = re.match(SIM_DIR_PATTERN, sim_dir)
	if match:
		timestamp = match.group(1)
		description = match.group(2).replace('_', ' ')
	else:
		timestamp = fp.timestamp()
		description = 'a manual run'

	return timestamp, description


class RunSimulation(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli simulation'

	def help(self):
		"""Return help text for the Command Line Interface."""
		return '''Run a {}.
				If the sim_path ends with a dir like
				"20190704.101500__Latest_sim_run", this will get the
				timestamp and description from the path to write into
				metadata.cPickle.
				The command line option names are long but you can use any
				unambiguous prefix.'''.format(self.description())

	def define_parameters(self, parser):
		def add_option(name, key, datatype, help):
			"""Add an option with the given name and datatype to the parser using
			DEFAULT_SIMULATION_KWARGS[key] for the default value.
			"""
			default = DEFAULT_SIMULATION_KWARGS[key]
			self.define_option(parser, name, datatype, default, help)
		def add_bool_option(name, key, help):
			"""Add a boolean option parameter with the given name to the parser
			using DEFAULT_SIMULATION_KWARGS[key] for the default value. The CLI
			input can be `--name`, `--no_name`, `--name true`, `--name false`,
			`--name 1`, `--name 0`, `--name=true`, etc.
			"""
			self.define_parameter_bool(
				parser, name, DEFAULT_SIMULATION_KWARGS[key], help)

		super(RunSimulation, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		parser.add_argument('-v', '--variant', nargs=3, default=['wildtype', '0', '0'],
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='The variant type name, first index, and last index. See'
				 ' models/ecoli/sim/variants/__init__.py for the possible'
				 ' variant choices. Default = wildtype 0 0')

		parser.add_argument('-g', '--generations', type=int, default=1,
			help='Number of cell generations to run. (Single daughters only.)'
				 ' Default = 1'
			)
		parser.add_argument('--total_gens', type=int,
			help='(int) Total number of generations to write into the'
				 ' metadata.cPickle file. Default = the value of --generations.')
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='First cell simulation seed. Default = 0'
			)
		self.define_option(parser, 'init_sims', int, 1,
			'Number of initial sims (seeds) per variant.')
		add_option('length_sec', 'lengthSec', int,
			help='The maximum simulation time, in seconds. Useful for short'
				 ' simulations; not so useful for multiple generations.'
				 ' Default is 3 hours'
			)
		add_option('timestep_safety_frac', 'timeStepSafetyFraction', float,
			help='Scale the time step by this factor if conditions are'
				 ' favorable, up the the limit of the max time step'
			)
		add_option('timestep_max', 'maxTimeStep', float,
			help='the maximum time step, in seconds'
			)
		add_option('timestep_update_freq', 'updateTimeStepFreq', int,
			help='frequency at which the time step is updated'  # TODO: explain
			)
		add_bool_option('mass_distribution', 'massDistribution',
			help='If true, a mass coefficient is drawn from a normal distribution'
				 ' centered on 1; otherwise it is set equal to 1'
			)
		add_bool_option('growth_rate_noise', 'growthRateNoise',
			help='If true, a growth rate coefficient is drawn from a normal'
				 ' distribution centered on 1; otherwise it is set equal to 1'
			)
		add_bool_option('d_period_division', 'dPeriodDivision',
			help='If true, ends simulation once D period has occurred after'
				 ' chromosome termination; otherwise simulation terminates once'
				 ' a given mass has been added to the cell'
			)
		add_bool_option('translation_supply', 'translationSupply',
			help='If true, the ribosome elongation rate is limited by the'
				 ' condition specific rate of amino acid supply; otherwise the'
				 ' elongation rate is set by condition')
		add_bool_option('variable_elongation_translation', 'variable_elongation_translation',
			help='if False, all elongation rates are the same for each type.'
				 'if True, elongation rates are faster for ribosomal proteins'
				 ' elongation rate is set by condition')
		add_bool_option('variable_elongation_transcription', 'variable_elongation_transcription',
			help='if False, all elongation rates are the same for each type.'
				 'if True, elongation rates are faster for RRNA transcripts'
				 ' elongation rate is set by condition')


	def parse_args(self):
		args = super(RunSimulation, self).parse_args()

		if args.total_gens is None:
			args.total_gens = args.generations

		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, 'simData_Fit_1.cPickle')
		fp.verify_file_exists(sim_data_file, 'Run runFitter?')

		timestamp, description = parse_timestamp_description(args.sim_path)

		variant_type = args.variant[0]
		variants_to_run = xrange(int(args.variant[1]), int(args.variant[2]) + 1)

		cli_sim_args = data.select_keys(vars(args), (
			'length_sec',
			'timestep_safety_frac',
			'timestep_max',
			'timestep_update_freq',
			'mass_distribution',
			'growth_rate_noise',
			'd_period_division',
			'translation_supply'))

		# Write the metadata file.

		cli_metadata_args = data.select_keys(vars(args), (
			'total_gens',
			'mass_distribution',
			'growth_rate_noise',
			'd_period_division',
			'translation_supply',
			'variable_elongation_translation',
			'variable_elongation_transcription'))

		metadata = dict(cli_metadata_args,
			git_hash=fp.run_cmdline("git rev-parse HEAD") or '--',
			git_branch=fp.run_cmdline("git symbolic-ref --short HEAD") or '--',
			description=description,
			time=timestamp,
			analysis_type=None,
			variant=variant_type,
			total_variants=str(len(variants_to_run)))

		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.SERIALIZED_METADATA_FILE)
		with open(metadata_path, "wb") as f:
			cPickle.dump(metadata, f, cPickle.HIGHEST_PROTOCOL)

		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i in variants_to_run:
			variant_directory = os.path.join(args.sim_path, variant_type + "_%06d" % i)
			variant_sim_data_directory = os.path.join(variant_directory, "kb")

			most_fit_filename = os.path.join(
				kb_directory, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME)
			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			fp.makedirs(variant_sim_data_directory)
			variant_metadata_directory = fp.makedirs(variant_directory,
				"metadata")
			task = VariantSimDataTask(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=most_fit_filename,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_directory,
				)
			task.run_task({})

			for j in xrange(args.seed, args.seed + args.init_sims):  # init sim seeds
				seed_directory = fp.makedirs(variant_directory, "%06d" % j)

				for k in xrange(args.generations):  # generation number k
					gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)

					# l is the daughter number among all of this generation's cells,
					# which is 0 for single-daughters but would span range(2**k) if
					# each parent had 2 daughters.
					l = 0
					cell_directory = fp.makedirs(gen_directory, "%06d" % l)
					cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

					options = dict(cli_sim_args,
						input_sim_data=variant_sim_data_modified_file,
						output_directory=cell_sim_out_directory,
						)

					if k == 0:
						task = SimulationTask(seed=j, **options)
					else:
						parent_gen_directory = os.path.join(seed_directory, "generation_%06d" % (k - 1))
						parent_cell_directory = os.path.join(parent_gen_directory, "%06d" % (l // 2))
						parent_cell_sim_out_directory = os.path.join(parent_cell_directory, "simOut")
						daughter_state_path = os.path.join(
							parent_cell_sim_out_directory, "Daughter%d" % (l % 2 + 1))
						task = SimulationDaughterTask(
							inherited_state_path=daughter_state_path,
							seed=(j + 1) * ((2**k - 1) + l),
							**options
							)
					task.run_task({})

if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
