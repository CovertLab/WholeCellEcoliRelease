"""
Common code for scripts that manually run analysis plots.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os
import sys

from wholecell.utils import constants, scriptBase, parallelization, filepath


class AnalysisBase(scriptBase.ScriptBase):
	"""Abstract base class for scripts that manually run analysis plots in an
	existing sim_dir: Defines a `sim_dir` command line parameter and sets some
	derived args.

	run() is still an abstract method.
	"""

	OUTPUT_SUBDIR = 'plotOut'  # recommended

	def __init__(self, analysis_plotter=None):
		"""Instantiate with an optional AnalysisPlot to run (this only uses its
		filename, not the instance), otherwise run all ACTIVE AnalysisPlots in
		the subclass' category.
		"""
		super(AnalysisBase, self).__init__()

		self.plot_name = None
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
					as defined in this category's __init__.py file. If omitted,
					the default is "CORE", which names the plots recommended
					for everyday development. Use "ACTIVE" to run all active
					plots in this category.''')

		parser.add_argument('-o', scriptBase.dashize('--output_prefix'),
			default='',
			help='Prefix for all the output plot filenames.')

		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. The given value will be'
				 ' limited to the available number of CPU cores. Default = 1.')

	def parse_args(self):
		"""Parse the command line args into an `argparse.Namespace`, including
		the `sim_dir` and `sim_path` args; sanitize args.plot; attach the
		`args.input_validation_data` path, the `args.metadata_path` path
		"<sim_path>/metadata/metadata.json", and the `args.metadata` dict
		loaded from `metadata_path`. If the superclass set `args.variant_dir`,
		also set `args.variant_dir_name` and metadata fields `variant_function`
		and `variant_index`.

		Overrides should first call super().
		"""
		args = super(AnalysisBase, self).parse_args()

		if self.plot_name:
			args.plot = [self.plot_name]

		args.input_validation_data = os.path.join(
			args.sim_path, 'kb', constants.SERIALIZED_VALIDATION_DATA)

		args.metadata_path = os.path.join(
			args.sim_path, 'metadata', constants.JSON_METADATA_FILE)
		args.metadata = filepath.read_json_file(args.metadata_path)

		if 'variant_dir' in args:
			args.variant_dir_name, variant_type, variant_index = args.variant_dir
			metadata = args.metadata
			metadata['variant_function'] = variant_type
			metadata['variant_index'] = variant_index

		args.cpus = parallelization.cpus(args.cpus)

		return args


class TestAnalysis(AnalysisBase):
	"""To test out the command line parser."""

	def __init__(self):
		super(TestAnalysis, self).__init__(analysis_plotter=self)

	def run(self, analysis_args):
		print("[TEST] Analysis args:", analysis_args)


if __name__ == '__main__':
	analysis = TestAnalysis()
	analysis.cli()
