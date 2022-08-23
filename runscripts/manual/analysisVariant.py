"""
Runs all variant analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import errno
import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisVariant import AnalysisVariantTask
from wholecell.utils import constants


class AnalysisVariant(AnalysisBase):
	"""Runs some or all the ACTIVE variant analysis plots for a given sim."""

	def define_parameters(self, parser):
		super().define_parameters(parser)
		self.define_path_selection(parser, 'variant', 'seed', 'generation')

	def update_args(self, args):
		super(AnalysisVariant, self).update_args(args)

		variant_dirs = self.list_variant_dirs(args.sim_path)  # list of tuples
		if not variant_dirs:
			raise IOError(errno.ENOENT,
				'No simulation variant directories found')

		metadata = args.metadata
		metadata['total_variants'] = str(len(variant_dirs))

	def run(self, args):
		output_dir = os.path.join(args.sim_path, constants.PLOTOUT_DIR)
		input_sim_data = os.path.join(args.sim_path,
			constants.KB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)

		task = AnalysisVariantTask(
			input_directory=args.sim_path,
			input_sim_data=input_sim_data,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix,
			variant_paths=args.variant_paths,
			seed_paths=args.seed_paths,
			generation_paths=args.generation_paths,
			only_successful=args.only_successful,
			**self.select_analysis_keys(args)
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisVariant()
	analysis.cli()
