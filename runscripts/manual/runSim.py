"""
Run a simple simulation, assuming you've run the Fitter. This does not run
multiple variants, multiple initial simulations, multiple generations, or
daughter simulations.

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
from wholecell.utils import constants, scriptBase
from wholecell.utils.filepath import makedirs


class RunSimulation(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def define_parameters(self, parser):
		super(RunSimulation, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		parser.add_argument('-v', '--variant', nargs=3,
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='The variant type name, first index, and last index. See'
				 ' models/ecoli/sim/variants/__init__.py for the possible'
				 ' variant choices. Default = wildtype 0 0')

	def parse_args(self):
		args = super(RunSimulation, self).parse_args()
		args.sim_path = scriptBase.find_sim_path(args.sim_dir)
		return args

	def run(self, args):
		sim_data_file = os.path.join(args.sim_path, 'kb', 'simData_Fit_1.cPickle')
		if not os.path.isfile(sim_data_file):
			raise IOError(errno.ENOENT,
				'Missing "{}".  Run the Fitter?'.format(sim_data_file))

		if args.variant:
			variant_type = args.variant[0]
			variants_to_run = range(args.variant[1], args.variant[2] + 1)
		else:
			variant_type = 'wildtype'
			variants_to_run = [0]

		print 'Variant {} {} .. {}'.format(
			variant_type, variants_to_run[0], variants_to_run[-1])

		# Set up variant, seed, and generation directories.
		# TODO(jerry): Skip most of this?
		for i in variants_to_run:
			variant_directory = makedirs(args.sim_path, variant_type + "_%06d" % i)
			variant_sim_data_directory = makedirs(variant_directory, "kb")
			variant_metadata_directory = makedirs(variant_directory, "metadata")
			variant_cohort_plot_directory = makedirs(variant_directory, "plotOut")

			j = 0  # init sim number
			seed_directory = makedirs(variant_directory, "%06d" % j)
			seed_plot_directory = makedirs(seed_directory, "plotOut")

			k = 0  # generation number
			gen_directory = makedirs(seed_directory, "generation_%06d" % k)

			l = 0  # daughter number
			cell_directory = makedirs(gen_directory, "%06d" % l)
			cell_sim_out_directory = makedirs(cell_directory, "simOut")
			cell_plot_out_directory = makedirs(cell_directory, "plotOut")

		# Write the metadata file.
		# TODO(jerry): Skip this? It's necessary but not sufficient for the
		# analysis scripts.
		metadata = {
			# "git_hash":           run_cmd(["git", "rev-parse", "HEAD"]),
			# "git_branch":         run_cmd(
			# 	["git", "symbolic-ref", "--short", "HEAD"]),
			# "git_diff":           run_cmd(["git", "diff"]),
			"description":        "a manual run",
			"time":               self.timestamp(),
			"total_gens":         1,
			"analysis_type":      None,
			"variant":            variant_type,
			# "mass_distribution":  MASS_DISTRIBUTION,
			# "growth_rate_noise":  GROWTH_RATE_NOISE,
			# "d_period_division":  D_PERIOD_DIVISION,
			# "translation_supply": TRANSLATION_SUPPLY,
		}
		metadata_dir = makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.SERIALIZED_METADATA_FILE)
		with open(metadata_path, "w") as f:
			cPickle.dump(metadata, f, cPickle.HIGHEST_PROTOCOL)

		# TODO(jerry): Optionally read/write to the
		# variant_directory/seed/generation/daughter/ path?
		task = SimulationTask(
			input_sim_data=sim_data_file,
			output_directory=makedirs(args.sim_path, 'sim'),
		)
		task.run_task({})


if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
