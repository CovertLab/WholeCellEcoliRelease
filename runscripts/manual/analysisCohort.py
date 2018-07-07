"""
Runs all cohort analysis plots for a given variant of a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisCohort import AnalysisCohortTask
from wholecell.utils import constants
from wholecell.utils import filepath


class AnalysisCohort(AnalysisBase):
	"""Runs some or all the ACTIVE cohort analysis plots for a given sim."""

	def define_parameters(self, parser):
		super(AnalysisCohort, self).define_parameters(parser)
		self.define_parameter_variant_index(parser)

	def parse_args(self):
		args = super(AnalysisCohort, self).parse_args()
		args.metadata['analysis_type'] = 'cohort'
		return args

	def run(self, args):
		sim_path = args.sim_path
		variant_dir_name = args.variant_dir_name

		input_variant_directory = os.path.join(sim_path, variant_dir_name)
		sim_data_modified = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		# TODO(jerry): Load simData_Modified into metadata?
		output_dir = filepath.makedirs(input_variant_directory, 'plotOut')

		task = AnalysisCohortTask(
			input_variant_directory=input_variant_directory,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			plots_to_run=args.plot,
			output_filename_prefix=args.output_prefix,
		)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisCohort()
	analysis.cli()
