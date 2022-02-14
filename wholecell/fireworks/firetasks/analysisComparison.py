"""
Run the plots for a given list of `models.ecoli.analysis.comparison` analyses,
or by default the ACTIVE plots listed in that package's `__init__.py`.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import annotations

from fireworks import explicit_serialize

from wholecell.fireworks.firetasks.analysisBase import AnalysisBase
import models.ecoli.analysis.comparison


@explicit_serialize
class AnalysisComparisonTask(AnalysisBase):

	_fw_name = "AnalysisComparisonTask"
	required_params = [
		"reference_sim_dir",
		"input_sim_dir",
		"output_plots_directory",
		"metadata",
		]
	optional_params = [
		"plot",  # absent or empty => run all active analysis plots
		"output_filename_prefix",
		"cpus",
		"compile",
		]
	MODULE_PATH = 'models.ecoli.analysis.comparison'
	TAGS = models.ecoli.analysis.comparison.TAGS

	def plotter_args(self, module_filename):
		self["metadata"] = dict(self["metadata"], analysis_type = "comparison")

		return (
			self["reference_sim_dir"],
			self["output_plots_directory"],
			self['output_filename_prefix'] + module_filename[:-3],
			self['input_sim_dir'],
			'',
			self["metadata"],
			)
