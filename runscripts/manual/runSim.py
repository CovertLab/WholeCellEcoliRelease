"""
Run a simple simulation, assuming you've run the Fitter. This does not run
multiple initial simulations, multiple generations, or multiple daughters per
generation.

TODO: Share lots of code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import errno
import os

from wholecell.fireworks.firetasks.simulation import SimulationTask
from wholecell.fireworks.firetasks import VariantSimDataTask
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

		parser.add_argument('-v', '--variant', nargs=3,
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='The variant type name, first index, and last index. See'
				 ' models/ecoli/sim/variants/__init__.py for the possible'
				 ' variant choices. Default = wildtype 0 0')

		add_option('length_sec', 'lengthSec', int,
			help='The maximum simulation time, in seconds. Useful for short'
				 ' simulations. Default is 3 hours'
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

	def parse_args(self):
		args = super(RunSimulation, self).parse_args()
		args.sim_path = scriptBase.find_sim_path(args.sim_dir)

		if not args.variant:
			args.variant = ['wildtype', '0', '0']

		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, 'simData_Fit_1.cPickle')
		if not os.path.isfile(sim_data_file):
			raise IOError(errno.ENOENT,
				'Missing "{}".  Run the Fitter?'.format(sim_data_file))

		variant_type = args.variant[0]
		variants_to_run = xrange(int(args.variant[1]), int(args.variant[2]) + 1)

		# Write the metadata file.
		metadata = {
			"git_hash":           fp.run_cmd(line="git rev-parse HEAD"),
			"git_branch":         fp.run_cmd(line="git symbolic-ref --short HEAD"),
			# "git_diff":           fp.run_cmd(line="git diff"),
			"description":        "a manual run",
			"time":               self.timestamp(),
			"total_gens":         1,
			"analysis_type":      None,
			"variant":            variant_type,
			"mass_distribution":  args.mass_distribution,
			"growth_rate_noise":  args.growth_rate_noise,
			"d_period_division":  args.d_period_division,
			"translation_supply": args.translation_supply,
			}
		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.SERIALIZED_METADATA_FILE)
		with open(metadata_path, "wb") as f:
			cPickle.dump(metadata, f, cPickle.HIGHEST_PROTOCOL)

		# TODO(jerry) Also write the redundant individual metadata_dir/key files?
		# If we want text metadata, just write metadata in Pickle protocol 0 or JSON.


		# Set up variant, seed, and generation directories.
		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i in variants_to_run:
			variant_directory = fp.makedirs(args.sim_path, variant_type + "_%06d" % i)
			variant_sim_data_directory = fp.makedirs(variant_directory, "kb")
			variant_metadata_directory = fp.makedirs(variant_directory, "metadata")

			most_fit_filename = os.path.join(
				kb_directory, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME)
			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			task = VariantSimDataTask(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=most_fit_filename,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_directory,
				)
			task.run_task({})

			j = 0  # init sim number. Don't loop over range(N_INIT_SIMS).
			seed_directory = fp.makedirs(variant_directory, "%06d" % j)

			k = 0  # generation number. (NOTE: Looping over range(N_GENS) would
			# require selecting SimulationTask vs. SimulationDaughterTask.)
			gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)

			l = 0  # daughter number. Don't support 2**k daughters.
			cell_directory = fp.makedirs(gen_directory, "%06d" % l)
			cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

			task = SimulationTask(
				input_sim_data=variant_sim_data_modified_file,
				output_directory=cell_sim_out_directory,
				seed=j,
				length_sec=args.length_sec,
				timestep_safety_frac=args.timestep_safety_frac,
				timestep_max=args.timestep_max,
				timestep_update_freq=args.timestep_update_freq,
				mass_distribution=args.mass_distribution,
				growth_rate_noise=args.growth_rate_noise,
				d_period_division=args.d_period_division,
				translation_supply=args.translation_supply,
				)
			task.run_task({})


if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
