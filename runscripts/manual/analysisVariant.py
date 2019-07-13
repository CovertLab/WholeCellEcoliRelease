"""
Runs all variant analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import errno
import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisVariant import AnalysisVariantTask


class AnalysisVariant(AnalysisBase):
	"""Runs some or all the ACTIVE variant analysis plots for a given sim."""

	def parse_args(self):
		args = super(AnalysisVariant, self).parse_args()

		variant_dirs = self.list_variant_dirs(args.sim_path)  # list of tuples
		if not variant_dirs:
			raise IOError(errno.ENOENT,
				'No simulation variant directories found')

		metadata = args.metadata
		metadata['total_variants'] = str(len(variant_dirs))

		return args

	def run(self, args):
		output_dir = os.path.join(args.sim_path, self.OUTPUT_SUBDIR)

		task = AnalysisVariantTask(
			input_directory=args.sim_path,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			plots_to_run=args.plot,
			output_filename_prefix=args.output_prefix,
			cpus=args.cpus,
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisVariant()
	analysis.cli()
