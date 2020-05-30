"""
Builds causality network for a given variant of a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.buildCausalityNetwork import BuildCausalityNetworkTask
from wholecell.utils import constants, filepath


CAUSALITY_ENV_VAR = 'CAUSALITY_SERVER'


class BuildCausalityNetwork(AnalysisBase):
	"""Builds causality network for a given sim."""

	def define_parameters(self, parser):
		super(BuildCausalityNetwork, self).define_parameters(parser)
		self.define_parameter_variant_index(parser)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='The initial simulation number (int). The value will get'
				 ' formatted as a subdirectory name like "000000". Default = 0.')
		parser.add_argument('-g', '--generation', type=int, default=0,
			help='The generation number (int). The value will get formatted'
				 ' as a subdirectory name like "generation_000000". Default = 0.')
		parser.add_argument('-d', '--daughter', type=int, default=0,
			help='The daughter number (int). The value will get formatted as'
				 ' a subdirectory name like "000000". Default = 0.')
		parser.add_argument('--check_sanity', action='store_true',
			help='Check network sanity.')
		parser.add_argument('-f', '--force', action='store_true',
			help='Forces a rebuild of the causality network if set.')
		parser.add_argument('--show', action='store_true',
			help='If set, attempts to show the causality visualization after'
			' processing data.')
		self.define_range_options(parser, 'variant', 'seed', 'generation')

	def update_args(self, args):
		super(BuildCausalityNetwork, self).update_args(args)

		args.seed_str = '%06d' % (args.seed,)
		args.gen_str = 'generation_%06d' % (args.generation,)
		args.daughter_str = '%06d' % (args.daughter,)

		metadata = args.metadata
		metadata['analysis_type'] = 'causality_network'
		metadata['seed'] = args.seed_str
		metadata['gen'] = args.gen_str

		return args

	def run(self, args):
		sim_path = args.sim_path
		variant_dir_name = args.variant_dir_name

		dirs = os.path.join(args.seed_str, args.gen_str, args.daughter_str)

		input_variant_directory = os.path.join(sim_path, variant_dir_name)
		input_dir = os.path.join(input_variant_directory, dirs, 'simOut')
		sim_data_modified = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		network_output_dir = os.path.join(input_variant_directory, 'kb')
		dynamics_output_dir = os.path.join(input_variant_directory, dirs, 'seriesOut')

		task = BuildCausalityNetworkTask(
			input_results_directory=input_dir,
			input_sim_data=sim_data_modified,
			output_network_directory=network_output_dir,
			output_dynamics_directory=dynamics_output_dir,
			check_sanity=args.check_sanity,
			metadata=args.metadata,
			force_update=args.force)
		task.run_task({})

		server_script = os.environ.get(CAUSALITY_ENV_VAR, '')
		if args.show and os.path.exists(server_script):
			# See #890 if running command fails due to differences in pyenv
			# versions - might need to cd to repo and activate pyenv
			cmd = 'python {} {} {}'.format(
				server_script, network_output_dir, dynamics_output_dir)
			print('Serving causality site with command:\n{}\n\nCtrl+c to exit.'
				.format(cmd))
			filepath.run_cmdline(cmd, timeout=None)
		else:
			print('\nNOTE: Use the --show flag and have ${0} set to your local'
				' path to site/server.py in the causality repo (eg. run'
				' `export {0}=~/path/to/causality/site/server.py` with your'
				' path) to automatically open the viewer for this data.\n'
				.format(CAUSALITY_ENV_VAR))


if __name__ == '__main__':
	analysis = BuildCausalityNetwork()
	analysis.cli()
