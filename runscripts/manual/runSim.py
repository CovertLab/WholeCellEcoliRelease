"""
Run a simple simulation, assuming you've run the parameter calculator first (Parca).
This does not run multiple initial simulations or multiple daughters per generation.

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import errno
import os

from wholecell.fireworks.firetasks import SimulationDaughterTask, SimulationTask, VariantSimDataTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp


class RunSimulation(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli simulation'

	def help(self):
		"""Return help text for the Command Line Interface."""
		return ('Run a {}. (The option names are long but you can use any'
				' unambiguous prefixes.)'.format(self.description()))

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
				 ' variant choices. Default = wildtype 0 0'
			)
		parser.add_argument('-g', '--generations', type=int, default=1,
			help='Number of cell generations to run. (Single daughters only.)'
				 ' Default = 1'
			)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='Cell simulation seed. Default = 0'
			)
		parser.add_argument('-t', '--timeline', type=str, default='0 minimal',
			help='set timeline. Default = "0 minimal". See'
				 ' environment/condition/make_media.py, make_timeline() for'
				 ' timeline formatting details'
			)
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
				 ' elongation rate is set by condition'
			)
		add_bool_option('trna_charging', 'trna_charging',
			help='if True, tRNA charging reactions are modeled and the ribosome'
				 ' elongation rate is set by the amount of charged tRNA	present.'
				 ' This option will override TRANSLATION_SUPPLY in the simulation.'
			)

	def parse_args(self):
		args = super(RunSimulation, self).parse_args()
		args.sim_path = scriptBase.find_sim_path(args.sim_dir)

		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		if not os.path.isfile(sim_data_file):
			raise IOError(errno.ENOENT,
				'Missing "{}".  Run the Parca?'.format(sim_data_file))

		variant_type = args.variant[0]
		variants_to_run = xrange(int(args.variant[1]), int(args.variant[2]) + 1)

		# Write the metadata file.
		metadata = {
			"git_hash":           fp.run_cmdline("git rev-parse HEAD") or '--',
			"git_branch":         fp.run_cmdline("git symbolic-ref --short HEAD") or '--',
			"description":        "a manual run",
			"time":               fp.timestamp(),
			"total_gens":         args.generations,
			"analysis_type":      None,
			"variant":            variant_type,
			"timeline":           args.timeline,
			"mass_distribution":  args.mass_distribution,
			"growth_rate_noise":  args.growth_rate_noise,
			"d_period_division":  args.d_period_division,
			"translation_supply": args.translation_supply,
			"trna_charging":      args.trna_charging,
			}
		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)


		# Set up variant, seed, and generation directories.
		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i in variants_to_run:
			variant_directory = fp.makedirs(args.sim_path, variant_type + "_%06d" % i)
			variant_sim_data_directory = fp.makedirs(variant_directory, "kb")
			variant_metadata_directory = fp.makedirs(variant_directory, "metadata")

			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			task = VariantSimDataTask(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=sim_data_file,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_directory,
				)
			task.run_task({})

			j = args.seed  # init sim number. This could loop over a range(N_INIT_SIMS).
			seed_directory = fp.makedirs(variant_directory, "%06d" % j)

			for k in xrange(args.generations):  # generation number k
				gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)

				# l is the daughter number among all of this generation's cells,
				# which is 0 for single-daughters but would span range(2**k) if
				# every parent had 2 daughters.
				l = 0
				cell_directory = fp.makedirs(gen_directory, "%06d" % l)
				cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

				options = dict(
					input_sim_data=variant_sim_data_modified_file,
					output_directory=cell_sim_out_directory,
					timeline=args.timeline,
					seed=j,
					length_sec=args.length_sec,
					timestep_safety_frac=args.timestep_safety_frac,
					timestep_max=args.timestep_max,
					timestep_update_freq=args.timestep_update_freq,
					mass_distribution=args.mass_distribution,
					growth_rate_noise=args.growth_rate_noise,
					d_period_division=args.d_period_division,
					translation_supply=args.translation_supply,
					trna_charging=args.trna_charging,
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
