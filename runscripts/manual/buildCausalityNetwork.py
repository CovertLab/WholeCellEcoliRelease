"""
Builds causality network for a given variant of a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os
import subprocess

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.buildCausalityNetwork import BuildCausalityNetworkTask
from wholecell.utils import constants


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
		sim_data_modified = os.path.join(input_variant_directory,
			constants.VKB_DIR,
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		dynamics_output_dir = os.path.join(input_variant_directory, dirs, 'seriesOut')

		task = BuildCausalityNetworkTask(
			input_results_directory=input_dir,
			input_sim_data=sim_data_modified,
			output_dynamics_directory=dynamics_output_dir,
			check_sanity=args.check_sanity,
			metadata=args.metadata)
		task.run_task({})

		# Optionally show the causality visualization.
		server_dir = os.environ.get(CAUSALITY_ENV_VAR, os.path.join('..', 'causality'))
		server_app = os.path.join('site', 'server.py')
		server_path = os.path.join(server_dir, server_app)
		if args.show and os.path.isfile(server_path):
			# See #890 if running command fails due to differences in pyenv
			# versions - might need to cd to repo and activate pyenv
			cmd = ['python', server_path, dynamics_output_dir]
			print(f'\nServing the Causality site via the command:\n  {cmd}\n'
				  f'Ctrl+C to exit.\n')
			subprocess.run(cmd)
		else:
			print('\nNOTE: Use the --show flag to automatically open the'
				' Casuality viewer on this data. You\'ll first need to'
				' `export {0}=~/path/to/causality` project unless the default'
				' (../causality) is good.\n'
				.format(CAUSALITY_ENV_VAR))


if __name__ == '__main__':
	analysis = BuildCausalityNetwork()
	analysis.cli()
