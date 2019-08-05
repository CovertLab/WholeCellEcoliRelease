"""
Run wcEcoli cell simulations, supporting multiple variants, multiple initial
seeds, and multiple generations, but only single daughters per generation.

Prerequisite: Run the parameter calculator (runParca.py).

Prerequisite: Generate the sim_data variant (makeVariants.py) before running
`runSim.py --require_variants`.

* Easy usage: runParca.py, then runSim.py, then analysis*.py.
* Fancy usage: runParca.py, makeVariants.py, lots of runSim.py and
  runDaughter.py runs, and analysis*.py, in a parallel workflow.

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

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
				metadata.json.
				The command line option names are long but you can use any
				unambiguous prefix.'''.format(self.description())

	def define_parameters(self, parser):

		super(RunSimulation, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		parser.add_argument(
			'-v', '--variant', nargs=3, default=['wildtype', '0', '0'],
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='''The variant type name, first index, and last index to make
				or require (depending on the --require_variants option). See
				models/ecoli/sim/variants/__init__.py for the variant
				type choices and their supported index ranges, e.g.: wildtype,
				condition, meneParams, metabolism_kinetic_objective_weight,
				nutrientTimeSeries, and param_sensitivity.
				Default = wildtype 0 0''')
		parser.add_argument(
			'-q', '--require-variants', type=bool, default=False,
			help='''true => require the sim_data variant(s) specified by the
				--variant option to already exist; false => make the variant(s).
				Run makeVariants.py to make sim_data variants.''')
		parser.add_argument(
			'-g', '--generations', type=int, default=1,
			help='Number of cell generations to run. (Single daughters only.)'
				 ' Default = 1')
		parser.add_argument(
			'-z', '--total-gens', type=int,
			help='(int) Total number of generations to write into the'
				 ' metadata.json file. Default = the value of --generations.')
		parser.add_argument(
			'-s', '--seed', type=int, default=0,
			help='First cell simulation seed. Default = 0')
		parser.add_argument(
			'-n', '--init-sims', type=int, default=1,
			help='Number of initial sims (seeds) per variant.')
		parser.add_argument(
			'-l', '--timeline', type=str, default='0 minimal',
			help='set timeline. Default = "0 minimal". See'
				 ' environment/condition/make_media.py, make_timeline() for'
				 ' timeline formatting details')
		parser.add_argument(
			'-t', '--length-sec', default=DEFAULT_SIMULATION_KWARGS['lengthSec'], type=int,
			help='The maximum simulation time, in seconds. Useful for short'
				 ' simulations; not so useful for multiple generations.'
				 ' Default is 3 hours')
		parser.add_argument(
			'-f', '--timestep-safety-frac',
			default=DEFAULT_SIMULATION_KWARGS['timeStepSafetyFraction'], type=float,
			help='Scale the time step by this factor if conditions are'
				 ' favorable, up the the limit of the max time step')
		parser.add_argument(
			'-x', '--timestep-max', default=DEFAULT_SIMULATION_KWARGS['maxTimeStep'], type=float,
			help='the maximum time step, in seconds')
		parser.add_argument(
			'-u', '--timestep-update-freq', default=DEFAULT_SIMULATION_KWARGS['updateTimeStepFreq'], type=int,
			help='frequency at which the time step is updated')  # TODO: explain
		parser.add_argument(
			'-m', '--mass-distribution', default=DEFAULT_SIMULATION_KWARGS['massDistribution'], type=bool,
			help='If true, a mass coefficient is drawn from a normal distribution'
				 ' centered on 1; otherwise it is set equal to 1')
		parser.add_argument(
			'-o', '--growth-rate-noise', default=DEFAULT_SIMULATION_KWARGS['growthRateNoise'], type=bool,
			help='If true, a growth rate coefficient is drawn from a normal'
				 ' distribution centered on 1; otherwise it is set equal to 1')
		parser.add_argument(
			'-d', '--d-period-division', default=DEFAULT_SIMULATION_KWARGS['dPeriodDivision'], type=bool,
			help='If true, ends simulation once D period has occurred after'
				 ' chromosome termination; otherwise simulation terminates once'
				 ' a given mass has been added to the cell')
		parser.add_argument(
			'-e', '--variable-elongation-transcription', default=DEFAULT_SIMULATION_KWARGS['variable_elongation_transcription'], type=bool,
			help='If true, runs various transcripts at different elongation rates')
		parser.add_argument(
			'-a', '--variable-elongation-translation', default=DEFAULT_SIMULATION_KWARGS['variable_elongation_translation'], type=bool,
			help='If true, translates various transcripts at different elongation rates')
		parser.add_argument(
			'-y', '--translation-supply', default=DEFAULT_SIMULATION_KWARGS['translationSupply'], type=bool,
			help='If true, the ribosome elongation rate is limited by the'
				 ' condition specific rate of amino acid supply; otherwise the'
				 ' elongation rate is set by condition')
		parser.add_argument(
			'-r', '--trna-charging', default=DEFAULT_SIMULATION_KWARGS['trna_charging'], type=bool,
			help='if True, tRNA charging reactions are modeled and the ribosome'
				 ' elongation rate is set by the amount of charged tRNA	present.'
				 ' This option will override TRANSLATION_SUPPLY in the simulation.')

	def parse_args(self):
		args = super(RunSimulation, self).parse_args()

		if args.total_gens is None:
			args.total_gens = args.generations

		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		fp.verify_file_exists(sim_data_file, 'Run runParca?')

		timestamp, description = parse_timestamp_description(args.sim_path)

		variant_type = args.variant[0]
		variant_spec = (variant_type, int(args.variant[1]), int(args.variant[2]))

		cli_sim_args = data.select_keys(vars(args),
			('timeline', 'length_sec', 'timestep_safety_frac', 'timestep_max',
			 'timestep_update_freq', 'mass_distribution', 'growth_rate_noise',
			 'd_period_division', 'translation_supply', 'trna_charging',
			 'variable_elongation_translation', 'variable_elongation_transcription'))

		# Write the metadata file.
		cli_metadata_args = data.select_keys(vars(args),
			('total_gens', 'timeline', 'mass_distribution', 'growth_rate_noise',
			'd_period_division', 'translation_supply', 'trna_charging'))
		metadata = dict(cli_metadata_args,
			git_hash=fp.run_cmdline("git rev-parse HEAD") or '--',
			git_branch=fp.run_cmdline("git symbolic-ref --short HEAD") or '--',
			description=description,
			time=timestamp,
			analysis_type=None,
			variant=variant_type,
			total_variants=str(variant_spec[2] + 1 - variant_spec[1]),
			)
		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)


		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i, subdir in fp.iter_variants(*variant_spec):
			variant_directory = os.path.join(args.sim_path, subdir)
			variant_sim_data_directory = os.path.join(variant_directory,
				VariantSimDataTask.OUTPUT_SUBDIR_KB)

			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			if args.require_variants:
				fp.verify_file_exists(
					variant_sim_data_modified_file, 'Run makeVariants?')
			else:
				variant_metadata_directory = os.path.join(variant_directory,
					VariantSimDataTask.OUTPUT_SUBDIR_METADATA)
				task = VariantSimDataTask(
					variant_function=variant_type,
					variant_index=i,
					input_sim_data=sim_data_file,
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
						seed=j,
						)

					if k == 0:
						task = SimulationTask(**options)
					else:
						parent_gen_directory = os.path.join(seed_directory, "generation_%06d" % (k - 1))
						parent_cell_directory = os.path.join(parent_gen_directory, "%06d" % (l // 2))
						parent_cell_sim_out_directory = os.path.join(parent_cell_directory, "simOut")
						daughter_state_path = os.path.join(
							parent_cell_sim_out_directory,
							constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
						task = SimulationDaughterTask(
							inherited_state_path=daughter_state_path,
							**options
							)
					task.run_task({})


if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
