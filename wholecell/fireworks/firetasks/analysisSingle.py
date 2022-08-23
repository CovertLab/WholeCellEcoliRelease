"""
Run the plots for a given list of `models.ecoli.analysis.single` analyses, or
by default the ACTIVE plots listed in that package's `__init__.py`.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import, division, print_function

from fireworks import explicit_serialize

from wholecell.fireworks.firetasks.analysisBase import AnalysisBase
import models.ecoli.analysis.single


@explicit_serialize
class AnalysisSingleTask(AnalysisBase):

	_fw_name = "AnalysisSingleTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]
	optional_params = [
		"plot",  # absent or empty => run all active analysis plots
		"output_filename_prefix",
		"cpus",
		"compile",
		]
	MODULE_PATH = 'models.ecoli.analysis.single'
	TAGS = models.ecoli.analysis.single.TAGS

	def plotter_args(self, module_filename):
		self["metadata"] = dict(self["metadata"], analysis_type = "single")

		return (
			self["input_results_directory"],
			self["output_plots_directory"],
			self['output_filename_prefix'] + module_filename[:-3],
			self["input_sim_data"],
			self["input_validation_data"],
			self["metadata"],
			)
