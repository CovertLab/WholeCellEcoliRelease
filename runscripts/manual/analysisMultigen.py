"""
Run all multigen analysis plots for initial sim #0 of a given variant of a
given simulation.

Run with '-h' for command line help.

Set the environment variable WC_ANALYZE_FAST to run multiple analysis scripts
in parallel processes.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisMultiGen import AnalysisMultiGenTask
from wholecell.utils import constants
from wholecell.utils import filepath


class AnalysisMultigen(AnalysisBase):
	"""Runs some or all the ACTIVE multigen analysis plots for a given sim."""

	def define_parameters(self, parser):
		super(AnalysisMultigen, self).define_parameters(parser)
		self.define_parameter_variant_index(parser)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='The initial simulation number (int). The value will get'
				 ' formatted as a subdirectory name like "000000". Default = 0.')

	def parse_args(self):
		args = super(AnalysisMultigen, self).parse_args()

		args.seed_str = '%06d' % (args.seed,)

		args.metadata['analysis_type'] = 'multigen'
		args.metadata["seed"] = args.seed_str
		return args

	def run(self, args):
		sim_path = args.sim_path
		variant_dir_name = args.variant_dir_name

		input_variant_directory = os.path.join(sim_path, variant_dir_name)
		input_path = os.path.join(input_variant_directory, args.seed_str)
		sim_data_modified = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		output_dir = filepath.makedirs(input_path, "plotOut")

		task = AnalysisMultiGenTask(
			input_seed_directory=input_path,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			plots_to_run=args.plot,
			output_filename_prefix=args.output_prefix,
		)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisMultigen()
	analysis.cli()
