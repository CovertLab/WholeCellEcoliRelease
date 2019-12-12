"""
Runs all cohort analysis plots for a given variant of a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisCohort import AnalysisCohortTask
from wholecell.utils import constants


class AnalysisCohort(AnalysisBase):
	"""Runs some or all the ACTIVE cohort analysis plots for a given sim."""

	def define_parameters(self, parser):
		super(AnalysisCohort, self).define_parameters(parser)
		self.define_parameter_variant_index(parser)
		self.define_range_options(parser, 'variant')

	def run(self, args):
		sim_path = args.sim_path
		variant_dir_name = args.variant_dir_name

		input_variant_directory = os.path.join(sim_path, variant_dir_name)
		sim_data_modified = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		output_dir = os.path.join(input_variant_directory, self.OUTPUT_SUBDIR)

		task = AnalysisCohortTask(
			input_variant_directory=input_variant_directory,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix,
			**self.select_analysis_keys(args)
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisCohort()
	analysis.cli()
