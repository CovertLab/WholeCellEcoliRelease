"""
Run a simple daughter simulation.  This does not run multiple generations.

Prerequisite: Run the parameter calculator (runParca.py).

Prerequisite: Run the parent generations (runSim.py, runDaughter.py).

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os

from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp


class RunDaughter(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli daughter simulation'

	def help(self):
		"""Return help text for the Command Line Interface."""
		return ('Run a {}. (The option names are long but you can use any'
				' unambiguous prefix.)'.format(self.description()))

	def define_parameters(self, parser):
		super(RunDaughter, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)
		self.define_parameter_variant_index(parser)

		parser.add_argument('-g', '--generation', type=int, default=1,
			help='Generation number to run. Default = 1')
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='Cell simulation seed. Default = 0')
		parser.add_argument('-d', '--daughter', type=int, default=0,
			help='Cell daughter number (only 0 if single_daughters). Default = 0')

		self.define_sim_options(parser)
		self.define_elongation_options(parser)

	def parse_args(self):
		args = super(RunDaughter, self).parse_args()
		args.sim_path = scriptBase.find_sim_path(args.sim_dir)

		return args

	def run(self, args):
		args.variant_dir_name, variant_type, variant_index = args.variant_dir
		input_variant_directory = os.path.join(args.sim_path, args.variant_dir_name)

		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		fp.verify_file_exists(sim_data_file, 'Run runParca?')

		cli_sim_args = data.select_keys(vars(args), scriptBase.SIM_KEYS)

		j = args.seed
		k = args.generation
		l = args.daughter

		# First generation not supported
		if k == 0:
			raise Exception('runDaughter is for daughter simulations. To run a first'
				' generation sim and for more options, use runSim')

		# Directory paths
		variant_sim_data_directory = fp.makedirs(input_variant_directory, "kb")
		variant_sim_data_modified_file = os.path.join(
			variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)
		seed_directory = fp.makedirs(input_variant_directory, "%06d" % j)
		gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)
		cell_directory = fp.makedirs(gen_directory, "%06d" % l)
		cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")
		parent_gen_directory = os.path.join(seed_directory, "generation_%06d" % (k - 1))
		parent_cell_directory = os.path.join(parent_gen_directory, "%06d" % (l // 2))
		parent_cell_sim_out_directory = os.path.join(parent_cell_directory, "simOut")
		daughter_state_path = os.path.join(
			parent_cell_sim_out_directory,
			constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))

		options = dict(cli_sim_args,
			input_sim_data=variant_sim_data_modified_file,
			output_directory=cell_sim_out_directory,
			seed=(j + 1) * ((2**k - 1) + l),
			)

		task = SimulationDaughterTask(
			inherited_state_path=daughter_state_path,
			**options
			)
		task.run_task({})


if __name__ == '__main__':
	script = RunDaughter()
	script.cli()
