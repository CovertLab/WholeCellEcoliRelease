"""
Common code for scripts that manually run analysis plots.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

import abc
import argparse
import os
import sys
from typing import Any, Dict, Optional

from models.ecoli.analysis.analysisPlot import AnalysisPlot
from wholecell.utils import constants, data, scriptBase, parallelization, filepath


class AnalysisBase(scriptBase.ScriptBase, metaclass=abc.ABCMeta):
	"""Abstract base class for scripts that manually run analysis plots in an
	existing sim_dir: Defines a `sim_dir` command line parameter and sets some
	derived args.

	run() is still an abstract method.
	"""

	def __init__(self, analysis_plotter=None):
		# type: (Optional[AnalysisPlot]) -> None
		"""Instantiate with an optional specific AnalysisPlot to run (this only
		uses its filename, not the instance), otherwise add the -p option to
		name the plots in the subclass' category to run, defaulting to the CORE
		list.
		"""
		super(AnalysisBase, self).__init__()

		self.plot_name = ''  # falsy => add_argument('-p')
		if analysis_plotter:
			module = sys.modules[type(analysis_plotter).__module__]
			self.plot_name = os.path.basename(module.__file__)

	def description(self):
		"""Describe the command line program."""
		return (self.plot_name + ' plot' if self.plot_name
				else type(self).__name__ + ' plots')

	def define_parameters(self, parser):
		"""Define command line parameters including `sim_dir` and `--plots`.

		When overriding, first call super().
		"""
		super(AnalysisBase, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		if not self.plot_name:  # not pre-specified
			parser.add_argument('-p', '--plot', nargs='+', default=[],
				help='''Names the analysis plots to run, e.g. plot filenames
					like "aaCounts.py" or "aaCounts" and tags like "METABOLISM"
					as defined in the __init__.py files. If omitted, the default
					is "DEFAULT", which will run the "CORE" tag with plots
					recommended for everyday development and the "VARIANT" tag
					for any variant specific plots. Use "ACTIVE" to run all
					active plots in this category. You can name specific
					analysis files but any analysis categories that don't have
					those files will print error messages.''')

		parser.add_argument('-o', scriptBase.dashize('--output_prefix'),
			default='',
			help='Prefix for all the output plot filenames.')

		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. The given value will be'
				 ' limited to the available number of CPU cores. Default = 1.')

		self.define_parameter_bool(parser, 'compile', False,
			'Compiles output images into one file (only for .png).')

	def select_analysis_keys(self, args):
		# type: (argparse.Namespace) -> Dict[str, Any]
		"""Select key/value pairs specific to analysis tasks"""
		return data.select_keys(vars(args), scriptBase.ANALYSIS_KEYS)

	def update_args(self, args):
		# type: (argparse.Namespace) -> None
		"""Update the command line args in an `argparse.Namespace`, including
		the `sim_dir` and `sim_path` args; sanitize args.plot; attach the
		`args.input_validation_data` path, the `args.metadata_path` path
		"<sim_path>/metadata/metadata.json", and the `args.metadata` dict
		loaded from `metadata_path`. If the superclass set `args.variant_dir`,
		also set `args.variant_dir_name` and metadata fields `variant_function`
		and `variant_index`.

		Overrides should first call super().
		"""
		super(AnalysisBase, self).update_args(args)

		if self.plot_name:
			args.plot = [self.plot_name]

		args.input_validation_data = os.path.join(
			args.sim_path, constants.KB_DIR, constants.SERIALIZED_VALIDATION_DATA)

		args.metadata_path = os.path.join(
			args.sim_path, constants.METADATA_DIR, constants.JSON_METADATA_FILE)
		args.metadata = filepath.read_json_file(args.metadata_path)

		if 'variant_dir' in args:
			args.variant_dir_name, variant_type, variant_index = args.variant_dir
			metadata = args.metadata
			metadata['variant_function'] = variant_type
			metadata['variant_index'] = variant_index

		args.cpus = parallelization.cpus(args.cpus)


class TestAnalysis(AnalysisBase):
	"""To test out the command line parser."""

	def __init__(self):
		self.plot = AnalysisPlot()
		super(TestAnalysis, self).__init__(analysis_plotter=self.plot)

	def run(self, analysis_args):
		print("[TEST] Analysis args:", analysis_args)


if __name__ == '__main__':
	analysis = TestAnalysis()
	analysis.cli()
