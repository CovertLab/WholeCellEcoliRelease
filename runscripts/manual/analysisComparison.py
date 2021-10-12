"""
Runs comparison analysis plots for two simulation output directories.

Run with '-h' for command line help.
"""

from __future__ import annotations

import errno
import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisComparison import AnalysisComparisonTask
from wholecell.utils import constants, scriptBase
from wholecell.utils import filepath as fp


class AnalysisComparison(AnalysisBase):
	"""Runs some or all the ACTIVE comparison analysis plots for a given sim."""

	def define_parameters(self, parser):
		super().define_parameters(parser)

		parser.add_argument('sim_dir2', nargs='?',
			help=f'''The reference simulation root directory to compare against.
				This argument can name a subdirectory of "out/" (optionally
				starting with "out/") or an absolute path. It defaults to
				sim_dir ± "{constants.OPERON_SUFFIX}".''')

	def update_args(self, args):
		super(AnalysisComparison, self).update_args(args)

		# args.sim_dir and args.sim_dir2 hold the CLI args; None == default.
		# args.sim_path and args.sim_path2 are the derived values that get used.
		if args.sim_dir2:
			args.sim_path2 = scriptBase.find_sim_path(args.sim_dir2)
		else:
			if args.sim_path.endswith(constants.OPERON_SUFFIX):
				args.sim_path2 = args.sim_path[:-len(constants.OPERON_SUFFIX)]
			else:
				args.sim_path2 = args.sim_path + constants.OPERON_SUFFIX
			fp.verify_dir_exists(args.sim_path2, "Failed to guess sim_dir2.")

		if args.sim_path == args.sim_path2:
			raise ValueError(
				f"SIM_PATH {args.sim_path} and SIM_PATH2 {args.sim_path2} must differ")

		variant_dirs = self.list_variant_dirs(args.sim_path)  # list of tuples
		if not variant_dirs:
			raise IOError(errno.ENOENT,
				'No simulation variant directories found')

		metadata = args.metadata
		metadata['total_variants'] = str(len(variant_dirs))

	def run(self, args):
		output_dir = os.path.join(args.sim_path, constants.COMPARISON_PLOTOUT_DIR)

		task = AnalysisComparisonTask(
			input_directory1=args.sim_path2,  # the reference
			input_directory2=args.sim_path,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix,
			**self.select_analysis_keys(args)
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisComparison()
	analysis.cli()
