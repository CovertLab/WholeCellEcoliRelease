"""
Runs all parca analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisParca import AnalysisParcaTask
from wholecell.utils import constants


class AnalysisParca(AnalysisBase):
	"""Runs some or all the ACTIVE parca analysis plots for a given sim."""

	def run(self, args):
		input_dir = os.path.join(args.sim_path, constants.KB_DIR)
		output_dir = os.path.join(args.sim_path, constants.KB_PLOT_OUTPUT_DIR)
		input_sim_data = os.path.join(input_dir, constants.SERIALIZED_SIM_DATA_FILENAME)

		task = AnalysisParcaTask(
			input_directory=input_dir,
			input_sim_data=input_sim_data,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix,
			**self.select_analysis_keys(args)
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisParca()
	analysis.cli()
